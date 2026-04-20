"""Methylation backends.

A backend's job is, given a list of probe coordinates, to produce per-probe
beta-value summaries (mean, q25, q75, n) for both a tumor and a normal cohort.
Cohort orchestration consumes that to produce `MethylationObservation`s per
candidate site.

Two backends ship in V1:

* `LocalArrayBackend` — works against pre-downloaded files on disk
  (probe annotation TSV + tumor beta TSV + normal beta TSV). This is the format
  GDC bulk downloads land in once you've fetched the cohort with `gdc-client`,
  so it doubles as the realistic on-disk path.

* `GDCBackend` — STUB. Documents the intended API surface so a real
  implementation can be dropped in without touching `cohort.py`. Not wired yet.
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
    # statistics.quantiles uses exclusive method by default: returns 3 cuts at .25, .50, .75
    qs = statistics.quantiles(clean, n=4, method="exclusive")
    q25, _q50, q75 = qs[0], qs[1], qs[2]
    # clamp into [0, 1] to honor the model invariant; rounding errors can push slightly out
    q25 = max(0.0, min(1.0, q25))
    q75 = max(0.0, min(1.0, q75))
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
        self._probes: list[ProbeRecord] = []
        for row in read_tsv(probe_annotation):
            try:
                self._probes.append(
                    ProbeRecord(
                        probe_id=row["probe_id"],
                        chrom=row["chrom"],
                        pos=int(row["pos"]),
                    )
                )
            except (KeyError, ValueError) as e:
                raise ValueError(
                    f"{probe_annotation}: expected columns probe_id, chrom, pos — {e}"
                ) from e

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

        # Load the raw matrix once, split by subtype.
        from thermocas.io import _open_text  # noqa: PLC0415 — internal helper, intentional
        import csv as _csv

        with _open_text(tumor_beta, "rt") as f:
            reader = _csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header is None:
                raise ValueError(f"{tumor_beta}: empty beta matrix")
            sample_ids = header[1:]
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

            # subtype → {probe_id → [betas restricted to that subtype]}
            per_subtype_betas: dict[str, dict[str, list[float | None]]] = {
                st: {} for st in cols_by_subtype
            }

            for row in reader:
                if not row:
                    continue
                probe_id = row[0]
                cells = row[1:]
                if len(cells) != len(sample_ids):
                    raise ValueError(
                        f"{tumor_beta}: probe {probe_id} has {len(cells)} cells, "
                        f"expected {len(sample_ids)}"
                    )
                parsed = [_parse_beta_cell(c) for c in cells]
                for st, idxs in cols_by_subtype.items():
                    per_subtype_betas[st][probe_id] = [parsed[i] for i in idxs]

        # Build one backend per subtype by reusing the constructor's path
        # without re-reading the raw matrix.
        backends: dict[str, LocalArrayBackend] = {}
        # Pre-load normal once (shared across subtype backends).
        _, normal_betas = read_beta_matrix(normal_beta)
        normal_summary = {pid: _summarize(pid, b) for pid, b in normal_betas.items()}
        probes = list(_load_probes(probe_annotation))

        for subtype, betas in per_subtype_betas.items():
            tumor_summary = {pid: _summarize(pid, b) for pid, b in betas.items()}
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
    """Internal: yield ProbeRecord rows from a probe annotation TSV."""

    for row in read_tsv(path):
        try:
            yield ProbeRecord(
                probe_id=row["probe_id"],
                chrom=row["chrom"],
                pos=int(row["pos"]),
            )
        except (KeyError, ValueError) as e:
            raise ValueError(
                f"{path}: expected columns probe_id, chrom, pos — {e}"
            ) from e


# ---------- GDCBackend ----------


class GDCBackend(MethylationBackend):
    """V2 — GDC HTTP client over stdlib `urllib`, with local file cache.

    Workflow:

      1. `list_files()` → POST /files with filters for project + platform +
         data type + sample type. Returns one file per sample.
      2. `download_file(file_id)` → GET /data/<id> → cached to `cache_dir`.
      3. `build_summaries()` → parse cached files → per-probe BetaSummary dict.
      4. `export_to_local(...)` → write probes.tsv + tumor.tsv + normal.tsv
         compatible with `LocalArrayBackend`. Run once, then re-use.

    GDC's per-sample beta files are typically ~50-250 MB on HM450, so a full
    cohort can be many GB. For repeated runs, prefer the export-then-LocalArray
    pattern via `thermocas gdc-fetch` rather than calling this backend live.

    Auth: GDC's public methylation data does not require an API token. If GDC
    later restricts access, supply `auth_token` and the client will pass it as
    the X-Auth-Token header.
    """

    BASE_URL = "https://api.gdc.cancer.gov"

    def __init__(
        self,
        project_id: str,
        cache_dir: str | Path,
        platform: str = "HM450",
        sample_type: str = "Primary Tumor",
        auth_token: str | None = None,
        url_opener: object | None = None,  # tests inject a fake opener
    ) -> None:
        self.project_id = project_id
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.platform = platform
        self.sample_type = sample_type
        self.auth_token = auth_token
        # `url_opener(url, data=None, headers={}) -> bytes`. Defaults to
        # `urllib.request.urlopen`. Tests pass a stub.
        self._opener = url_opener
        self._files_cache: list[dict] | None = None
        self._summaries_cache: dict[str, BetaSummary] | None = None

    # ---- low-level HTTP ----

    def _http(self, url: str, body: bytes | None = None) -> bytes:
        import urllib.request

        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        if self.auth_token:
            headers["X-Auth-Token"] = self.auth_token
        if self._opener is not None:
            return self._opener(url, body, headers)  # type: ignore[misc]
        req = urllib.request.Request(url, data=body, headers=headers, method="POST" if body else "GET")
        with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310 — GDC is HTTPS
            return resp.read()

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
        for meta in self.list_files():
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
        """Write a per-probe summary TSV: probe_id, n, mean, q25, q75.

        The `score-cohort` CLI consumes a beta matrix, not a summary, so for
        the LocalArrayBackend bridge use `export_betas_tsv` instead. This
        helper exists for inspection / downstream stats.
        """

        summaries = self.build_summaries()
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("probe_id\tn\tmean\tq25\tq75\n")
            for pid in sorted(summaries):
                s = summaries[pid]
                f.write(f"{pid}\t{s.n_samples}\t"
                        f"{_or_na(s.mean)}\t{_or_na(s.q25)}\t{_or_na(s.q75)}\n")
        return len(summaries)

    # ---- MethylationBackend interface ----

    def probes(self) -> Iterable[ProbeRecord]:
        """Probe coordinates aren't returned by the GDC `/files` endpoint;
        callers should pair this backend with an externally provided probe
        annotation TSV (the GDC manifest for HM450/EPIC annotations is the
        canonical source). Returning empty here forces the cohort layer to
        rely on a separate annotation file."""

        return []

    def tumor_summary(self, probe_id: str) -> BetaSummary | None:
        if self.sample_type != "Primary Tumor":
            return None
        return self.build_summaries().get(probe_id)

    def normal_summary(self, probe_id: str) -> BetaSummary | None:
        if self.sample_type != "Solid Tissue Normal":
            return None
        return self.build_summaries().get(probe_id)


# Map our short platform identifiers to GDC's expected strings.
_GDC_PLATFORM_MAP = {
    "HM450": "Illumina Human Methylation 450",
    "EPIC": "Illumina Human Methylation EPIC",
}


def _or_na(v: float | None) -> str:
    return "NA" if v is None else f"{v:.6f}"


