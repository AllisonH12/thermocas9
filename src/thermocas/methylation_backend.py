"""Methylation backends.

A `MethylationBackend` implementation, given a list of probe coordinates,
produces per-probe beta-value summaries (mean, q25, q75, n) for both a tumor
and a normal cohort. Cohort orchestration consumes that to produce
`MethylationObservation`s per candidate site.

Two `MethylationBackend` implementations ship:

* `LocalArrayBackend` — works against pre-downloaded raw beta matrices on
  disk (probe annotation TSV + tumor beta TSV + normal beta TSV, one column
  per sample). Use when you already have the per-sample matrix.

* `LocalSummaryBackend` — works against per-probe summary TSVs
  (probe_id, n, mean, q25, q75) plus a probe annotation TSV. This is the
  format `gdc-fetch` writes directly, so it's the canonical live-data path:
  `thermocas gdc-fetch → thermocas score-cohort --backend summary`.

One fetcher (not a `MethylationBackend`) ships:

* `GDCBackend` — stdlib-urllib client against the live GDC API with on-disk
  caching. Lists files, downloads them, aggregates per-probe, exports a
  summary TSV. Deliberately *not* a `MethylationBackend` subclass: live GDC
  traffic can't satisfy the backend contract because per-sample files don't
  carry probe coordinates and one fetcher only handles one sample-type side.
  Use it to materialize a cohort, then point `LocalSummaryBackend` at the
  output directory.
"""

from __future__ import annotations

import json
import statistics
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from thermocas.evidence import ProbeRecord
from thermocas.io import read_beta_matrix, read_sample_subtypes, read_tsv


@dataclass(frozen=True)
class BetaSummary:
    """Per-probe beta-value summary across one cohort's samples."""

    probe_id: str
    n_samples: int
    mean: float | None
    q25: float | None
    q75: float | None


class MethylationBackend(ABC):
    """Abstract methylation backend.

    A backend exposes:
      * `probes()`               → all available probes (with chrom + position)
      * `tumor_summary(probe)`   → BetaSummary in the tumor cohort
      * `normal_summary(probe)`  → BetaSummary in the normal cohort
    """

    @abstractmethod
    def probes(self) -> Iterable[ProbeRecord]:
        """Yield every assayed probe with chromosome and forward-strand position."""

    @abstractmethod
    def tumor_summary(self, probe_id: str) -> BetaSummary | None:
        """Return tumor-side BetaSummary for a probe, or None if missing."""

    @abstractmethod
    def normal_summary(self, probe_id: str) -> BetaSummary | None:
        """Return normal-side BetaSummary for a probe, or None if missing."""


# ---------- LocalArrayBackend ----------


def _summarize(probe_id: str, betas: list[float | None]) -> BetaSummary:
    """Compute mean + q25/q75 from a per-probe vector of beta values.

    NaN / missing samples are dropped before summarization. Returns None
    statistics when fewer than 2 valid samples (quantiles undefined).
    """

    clean = [b for b in betas if b is not None]
    n = len(clean)
    if n == 0:
        return BetaSummary(probe_id=probe_id, n_samples=0, mean=None, q25=None, q75=None)
    mean = statistics.fmean(clean)
    if n < 2:
        # quantiles undefined for n < 2; fall back to point summary
        return BetaSummary(probe_id=probe_id, n_samples=n, mean=mean, q25=mean, q75=mean)
    # `inclusive` (rather than exclusive) keeps quartiles within the empirical
    # support — important for small subtype splits where exclusive quartile
    # extrapolation invents methylation outside any observed sample.
    qs = statistics.quantiles(clean, n=4, method="inclusive")
    q25, _q50, q75 = qs[0], qs[1], qs[2]
    # Hard clamp to the observed sample range — defensive, since `inclusive`
    # already keeps us inside but we want the invariant to hold even if the
    # underlying stdlib changes its semantics.
    obs_min, obs_max = min(clean), max(clean)
    q25 = max(obs_min, min(obs_max, q25))
    q75 = max(obs_min, min(obs_max, q75))
    mean = max(q25, min(q75, mean))  # ensure mean lies within [q25, q75]
    return BetaSummary(probe_id=probe_id, n_samples=n, mean=mean, q25=q25, q75=q75)


class LocalArrayBackend(MethylationBackend):
    """Methylation backend over local TSV files.

    Expected file formats:

    Probe annotation TSV (header required):

        probe_id    chrom    pos
        cg00000029  chr16    53434200
        cg00000108  chr3     37459206

    Beta matrix TSV (header required, first column is probe ID):

        probe_id    sample1   sample2   sample3
        cg00000029  0.85      0.91      0.88
        cg00000108  0.45      NA        0.42

    NA / blank / NaN cells are treated as missing.
    """

    def __init__(
        self,
        probe_annotation: str | Path,
        tumor_beta: str | Path,
        normal_beta: str | Path,
    ) -> None:
        # `_load_probes` enforces the duplicate-probe_id invariant + rejects
        # malformed rows; routing both LocalArrayBackend and LocalSummaryBackend
        # through it keeps the two backends' annotation contract identical.
        self._probes = list(_load_probes(probe_annotation))

        _, tumor_betas = read_beta_matrix(tumor_beta)
        _, normal_betas = read_beta_matrix(normal_beta)

        self._tumor_summary = {pid: _summarize(pid, b) for pid, b in tumor_betas.items()}
        self._normal_summary = {pid: _summarize(pid, b) for pid, b in normal_betas.items()}

    def probes(self) -> Iterable[ProbeRecord]:
        return list(self._probes)

    def tumor_summary(self, probe_id: str) -> BetaSummary | None:
        return self._tumor_summary.get(probe_id)

    def normal_summary(self, probe_id: str) -> BetaSummary | None:
        return self._normal_summary.get(probe_id)

    @classmethod
    def split_by_subtype(
        cls,
        probe_annotation: str | Path,
        tumor_beta: str | Path,
        normal_beta: str | Path,
        sample_subtypes: str | Path,
    ) -> dict[str, "LocalArrayBackend"]:
        """V2 — return one `LocalArrayBackend` per subtype.

        Splits the tumor matrix on subtype labels and returns a mapping
        `subtype → backend`. The normal matrix is shared across all backends
        because in TCGA the normal compartment is rarely subtyped (the subtypes
        belong to the cancer side); callers can pre-split normals and pass per-
        subtype backends manually if needed.

        Subtype YAML schema:

            sample_id    subtype
            TCGA-AA-1234 LumA
            TCGA-BB-5678 LumB

        Use the resulting backends with `cohort_name` distinguished per
        subtype, e.g. by writing each output as `scored.BRCA-LumA.jsonl`. The
        pan-cancer aggregator treats each as its own cohort.
        """

        subtype_map = read_sample_subtypes(sample_subtypes)
        if not subtype_map:
            raise ValueError(f"{sample_subtypes}: no sample → subtype rows found")

        # Reuse the hardened `read_beta_matrix` loader instead of manually
        # reparsing the TSV. This inherits duplicate-probe-row rejection,
        # duplicate-sample-ID rejection, and cell-validity checks — the
        # V2 manual-parse path missed all three.
        sample_ids, tumor_betas = read_beta_matrix(tumor_beta)

        unknown = [s for s in sample_ids if s not in subtype_map]
        if unknown:
            raise ValueError(
                f"{tumor_beta}: {len(unknown)} sample(s) missing from "
                f"subtype map (e.g. {unknown[:3]})"
            )

        # subtype → list[col_idx] within the original matrix
        cols_by_subtype: dict[str, list[int]] = {}
        for i, sid in enumerate(sample_ids):
            cols_by_subtype.setdefault(subtype_map[sid], []).append(i)

        # Normal side shared across subtype backends (V2 design).
        _, normal_betas = read_beta_matrix(normal_beta)
        normal_summary = {pid: _summarize(pid, b) for pid, b in normal_betas.items()}
        probes = list(_load_probes(probe_annotation))

        backends: dict[str, LocalArrayBackend] = {}
        for subtype, idxs in cols_by_subtype.items():
            tumor_summary = {
                pid: _summarize(pid, [betas[i] for i in idxs])
                for pid, betas in tumor_betas.items()
            }
            backends[subtype] = cls.__new__(cls)
            backends[subtype]._probes = probes
            backends[subtype]._tumor_summary = tumor_summary
            backends[subtype]._normal_summary = normal_summary

        return backends


def _parse_beta_cell(cell: str) -> float | None:
    """Internal: parse a single beta cell as float in [0,1] or None for missing."""

    cell = cell.strip()
    if cell == "" or cell.upper() in {"NA", "NAN", "NULL"}:
        return None
    v = float(cell)
    if v < 0.0 or v > 1.0:
        raise ValueError(f"beta value {v} not in [0, 1]")
    return v


def _load_probes(path: str | Path) -> Iterable[ProbeRecord]:
    """Public-ish helper: yield ProbeRecord rows from a probe annotation TSV.

    Rejects duplicate probe_ids. Tumor/normal summaries are keyed only by
    probe_id; allowing the same probe_id at two different coordinates would
    alias one assay across multiple loci (fabricating EXACT evidence at
    whichever coordinate the classifier happens to pick).

    Tolerated by both LocalArrayBackend and LocalSummaryBackend.
    """

    seen: set[str] = set()
    for row in read_tsv(path):
        try:
            probe_id = row["probe_id"]
            if probe_id in seen:
                raise ValueError(
                    f"{path}: duplicate probe_id {probe_id!r} — annotations "
                    "must map each probe to exactly one genomic position"
                )
            seen.add(probe_id)
            yield ProbeRecord(
                probe_id=probe_id,
                chrom=row["chrom"],
                pos=int(row["pos"]),
            )
        except (KeyError, ValueError) as e:
            raise ValueError(
                f"{path}: expected columns probe_id, chrom, pos — {e}"
            ) from e


# ---------- LocalSummaryBackend ----------


class LocalSummaryBackend(MethylationBackend):
    """V3 — backend over per-probe summary TSVs (the format `gdc-fetch` writes).

    Inputs:
      * `probe_annotation` — TSV with probe_id, chrom, pos columns
      * `tumor_summary` and `normal_summary` — TSVs with
        `probe_id, n, mean, q25, q75` (NA-tolerant)

    The summary format is what `GDCBackend.export_summaries()` produces, so the
    intended workflow is:

        thermocas gdc-fetch --project TCGA-BRCA --output-dir cohort/
        thermocas score-cohort --backend summary \\
            --probe-annotation cohort/probes.tsv \\
            --tumor-summary cohort/tumor_summary.tsv \\
            --normal-summary cohort/normal_summary.tsv \\
            ...

    The class implements `MethylationBackend`, unlike `GDCBackend` which is
    just a fetcher.
    """

    def __init__(
        self,
        probe_annotation: str | Path,
        tumor_summary: str | Path,
        normal_summary: str | Path,
    ) -> None:
        self._probes = list(_load_probes(probe_annotation))
        self._tumor_summary = _load_summary_tsv(tumor_summary)
        self._normal_summary = _load_summary_tsv(normal_summary)

    def probes(self) -> Iterable[ProbeRecord]:
        return list(self._probes)

    def tumor_summary(self, probe_id: str) -> BetaSummary | None:
        return self._tumor_summary.get(probe_id)

    def normal_summary(self, probe_id: str) -> BetaSummary | None:
        return self._normal_summary.get(probe_id)


def _load_summary_tsv(path: str | Path) -> dict[str, BetaSummary]:
    """Parse a per-probe summary TSV (probe_id, n, mean, q25, q75 — NA-tolerant)."""

    out: dict[str, BetaSummary] = {}
    for row in read_tsv(path):
        try:
            n = int(row["n"])
        except (KeyError, ValueError) as e:
            raise ValueError(
                f"{path}: expected columns probe_id, n, mean, q25, q75 — {e}"
            ) from e
        probe_id = row["probe_id"]
        if probe_id in out:
            # See read_beta_matrix: probe_id is the join key; a duplicate row
            # is malformed input that would silently rewrite the earlier
            # summary. Reject rather than overwrite.
            raise ValueError(
                f"{path}: duplicate probe_id {probe_id!r} (last-wins overwrite "
                "would silently change scoring)"
            )
        out[probe_id] = BetaSummary(
            probe_id=probe_id,
            n_samples=n,
            mean=_parse_summary_cell(row.get("mean", "")),
            q25=_parse_summary_cell(row.get("q25", "")),
            q75=_parse_summary_cell(row.get("q75", "")),
        )
    return out


def _parse_summary_cell(cell: str) -> float | None:
    """Like _parse_beta_cell but accepts the rounded values export_summaries writes."""

    cell = cell.strip()
    if cell == "" or cell.upper() in {"NA", "NAN", "NULL"}:
        return None
    v = float(cell)
    # Be lenient with floating-point overshoot (export_summaries rounds to 6dp).
    if -1e-9 <= v <= 1.0 + 1e-9:
        return max(0.0, min(1.0, v))
    raise ValueError(f"summary value {v} not in [0, 1]")


# ---------- GDCBackend (fetcher only) ----------


class GDCBackend:
    """V2 — GDC HTTP fetcher over stdlib `urllib`, with on-disk file cache.

    **Note: this class is a fetcher, not a `MethylationBackend`.** It is
    intentionally not a `MethylationBackend` subclass: live GDC traffic
    cannot satisfy the backend contract because (a) the per-sample beta-value
    files don't carry probe coordinates, and (b) one instance only handles
    one sample-type side. The supported live-data path is:

        thermocas gdc-fetch --project TCGA-BRCA --output-dir cohort/
        # then:
        thermocas score-cohort --backend summary --probe-annotation cohort/probes.tsv \\
            --tumor-summary cohort/tumor_summary.tsv \\
            --normal-summary cohort/normal_summary.tsv \\
            ...

    The fetcher exposes:
      1. `list_files()` → POST /files with project + platform + data type +
         sample type filters; returns one record per sample.
      2. `download_file(file_id)` → GET /data/<id>, cached to `cache_dir`.
      3. `build_summaries()` → parse cached files into a per-probe BetaSummary
         dict (handles NA, capped by `max_files` for testing).
      4. `export_summaries(output_path)` → write a `probe_id, n, mean, q25, q75`
         TSV that `LocalSummaryBackend` consumes.

    Auth: GDC's public methylation data does not require an API token. Supply
    `auth_token` if a future endpoint restricts access; it'll be sent as
    `X-Auth-Token`.
    """

    BASE_URL = "https://api.gdc.cancer.gov"

    def __init__(
        self,
        project_id: str,
        cache_dir: str | Path,
        platform: str = "HM450",
        sample_type: str = "Primary Tumor",
        auth_token: str | None = None,
        max_files: int | None = None,
        url_opener: object | None = None,  # tests inject a fake opener
    ) -> None:
        self.project_id = project_id
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.platform = platform
        self.sample_type = sample_type
        self.auth_token = auth_token
        self.max_files = max_files  # cap files for testing or quota control
        # `url_opener(url, data=None, headers={}) -> bytes`. Defaults to
        # `urllib.request.urlopen`. Tests pass a stub.
        self._opener = url_opener
        self._files_cache: list[dict] | None = None
        self._summaries_cache: dict[str, BetaSummary] | None = None

    # ---- low-level HTTP ----

    def _http(self, url: str, body: bytes | None = None) -> bytes:
        """One HTTP call with retries for transient failures.

        GDC's /data endpoint occasionally times out on individual file
        downloads — bare urllib.urlopen with a 60s timeout was killing
        long-running gdc-fetch jobs after a single transient hiccup. We retry
        up to 3 times with exponential backoff + bump the per-call timeout
        to 300s so streaming a ~13 MB file over a slow link succeeds.
        """

        import time
        import urllib.request

        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        if self.auth_token:
            headers["X-Auth-Token"] = self.auth_token
        if self._opener is not None:
            return self._opener(url, body, headers)  # type: ignore[misc]

        last_err: Exception | None = None
        for attempt in range(3):
            try:
                req = urllib.request.Request(
                    url, data=body, headers=headers,
                    method="POST" if body else "GET",
                )
                with urllib.request.urlopen(req, timeout=300) as resp:  # noqa: S310 — GDC is HTTPS
                    return resp.read()
            except (TimeoutError, OSError) as e:
                last_err = e
                if attempt < 2:
                    time.sleep(2 ** attempt)  # 1s, then 2s
                    continue
                raise
        raise RuntimeError("unreachable") from last_err

    # ---- file listing ----

    def list_files(self) -> list[dict]:
        """Query GDC `/files` for one cohort's methylation files.

        Returns a list of file metadata dicts: at minimum `id`, `submitter_id`,
        `cases.0.case_id`. Cached in-memory after first call.
        """

        if self._files_cache is not None:
            return self._files_cache

        filters = {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": [self.project_id]}},
                {"op": "=", "content": {"field": "data_type", "value": ["Methylation Beta Value"]}},
                {"op": "=", "content": {"field": "platform", "value": [_GDC_PLATFORM_MAP.get(self.platform, self.platform)]}},
                {"op": "=", "content": {"field": "cases.samples.sample_type", "value": [self.sample_type]}},
            ],
        }
        body = json.dumps({
            "filters": filters,
            "fields": "file_id,file_name,cases.case_id,cases.samples.sample_type",
            "format": "JSON",
            "size": 10_000,
        }).encode()

        raw = self._http(f"{self.BASE_URL}/files", body=body)
        data = json.loads(raw)
        hits = data.get("data", {}).get("hits", [])
        self._files_cache = hits
        return hits

    # ---- per-file download (cached) ----

    def download_file(self, file_id: str) -> Path:
        """Download a single beta-value file to the cache directory; return path."""

        cached = self.cache_dir / f"{file_id}.txt"
        if cached.exists() and cached.stat().st_size > 0:
            return cached
        raw = self._http(f"{self.BASE_URL}/data/{file_id}")
        cached.write_bytes(raw)
        return cached

    # ---- per-probe summaries ----

    def build_summaries(self) -> dict[str, BetaSummary]:
        """Download every cohort file and aggregate into per-probe BetaSummary.

        File format (GDC harmonized methylation): TSV with no header, columns
        roughly `Composite Element REF\\tBeta_value\\t...`. We use only the
        first two columns. Missing values appear as `NA` or empty.
        """

        if self._summaries_cache is not None:
            return self._summaries_cache

        per_probe: dict[str, list[float | None]] = {}
        files = self.list_files()
        if self.max_files is not None:
            files = files[: self.max_files]
        for meta in files:
            file_id = meta.get("file_id") or meta.get("id")
            if not file_id:
                continue
            path = self.download_file(file_id)
            for row in path.read_text().splitlines():
                if not row or row.startswith("#"):
                    continue
                parts = row.split("\t")
                if len(parts) < 2:
                    continue
                probe_id = parts[0]
                if probe_id == "Composite Element REF":  # header row, skip
                    continue
                cell = parts[1].strip()
                if cell == "" or cell.upper() in {"NA", "NAN", "NULL"}:
                    per_probe.setdefault(probe_id, []).append(None)
                else:
                    try:
                        v = float(cell)
                    except ValueError:
                        continue
                    if 0.0 <= v <= 1.0:
                        per_probe.setdefault(probe_id, []).append(v)

        self._summaries_cache = {pid: _summarize(pid, betas) for pid, betas in per_probe.items()}
        return self._summaries_cache

    # ---- export to LocalArrayBackend-compatible TSVs ----

    def export_summaries(self, output_path: str | Path) -> int:
        """Write a per-probe summary TSV: `probe_id, n, mean, q25, q75`.

        This is the canonical post-fetch artifact — `LocalSummaryBackend`
        consumes it directly. Pair with a probe annotation TSV (chrom + pos)
        to drive `score-cohort --backend summary`.
        """

        summaries = self.build_summaries()
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("probe_id\tn\tmean\tq25\tq75\n")
            for pid in sorted(summaries):
                s = summaries[pid]
                f.write(f"{pid}\t{s.n_samples}\t"
                        f"{_or_na(s.mean)}\t{_or_na(s.q25)}\t{_or_na(s.q75)}\n")
        return len(summaries)


# Map our short platform identifiers to GDC's expected strings.
_GDC_PLATFORM_MAP = {
    "HM450": "Illumina Human Methylation 450",
    "EPIC": "Illumina Human Methylation EPIC",
}


def _or_na(v: float | None) -> str:
    return "NA" if v is None else f"{v:.6f}"


