"""IO helpers — JSONL, FASTA, TSV. Stdlib only.

The framework persists records as JSONL (one Pydantic model per line) so
catalogs, scored cohorts, and pan-cancer aggregates can stream and concatenate
without loading entire results into memory.
"""

from __future__ import annotations

import csv
import gzip
import io as _io
import json
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import IO, TypeVar

from pydantic import BaseModel

T = TypeVar("T", bound=BaseModel)


# ---------- JSONL ----------


def write_jsonl(path: str | Path, records: Iterable[BaseModel]) -> int:
    """Write Pydantic models to JSONL. Returns number of records written."""

    n = 0
    with _open_text(path, "wt") as f:
        for r in records:
            f.write(r.model_dump_json())
            f.write("\n")
            n += 1
    return n


def read_jsonl(path: str | Path, model: type[T]) -> Iterator[T]:
    """Stream Pydantic models from a JSONL file (gzip-aware)."""

    with _open_text(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield model.model_validate_json(line)


# ---------- FASTA ----------


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Stream (chrom, sequence) pairs from a FASTA file. Gzip-aware.

    Concatenates wrapped sequence lines per record. Strips description after
    the chromosome name (`>chr1 some description` → `chr1`).
    """

    with _open_text(path, "rt") as f:
        chrom: str | None = None
        chunks: list[str] = []
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(chunks)
                # Take the first whitespace-delimited token as the chromosome name.
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if chrom is not None:
            yield chrom, "".join(chunks)


# ---------- TSV ----------


def read_tsv(path: str | Path) -> Iterator[dict[str, str]]:
    """Stream rows from a tab-separated file with a header. Gzip-aware."""

    with _open_text(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        yield from reader


def read_beta_matrix(path: str | Path) -> tuple[list[str], dict[str, list[float | None]]]:
    """Load a probe × sample beta matrix.

    First column header is the probe-ID column (any name); remaining columns
    are sample IDs. Cell values are floats in [0, 1] or empty / 'NA' / 'NaN'
    for missing.

    Returns (sample_ids, {probe_id → [beta per sample]}).
    """

    samples: list[str] = []
    probe_betas: dict[str, list[float | None]] = {}

    with _open_text(path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{path}: empty beta matrix")
        samples = header[1:]
        for row in reader:
            if not row:
                continue
            probe_id = row[0]
            if probe_id in probe_betas:
                # probe_id is the join key into annotation + scoring; a duplicate
                # row is malformed input. Reject rather than silently overwriting
                # the earlier vector.
                raise ValueError(
                    f"{path}: duplicate probe_id {probe_id!r} (last-wins overwrite "
                    "would silently change cohort summaries)"
                )
            betas: list[float | None] = []
            for cell in row[1:]:
                cell = cell.strip()
                if cell == "" or cell.upper() in {"NA", "NAN", "NULL"}:
                    betas.append(None)
                    continue
                v = float(cell)
                if v < 0.0 or v > 1.0:
                    raise ValueError(
                        f"{path}: probe {probe_id} beta {v} not in [0, 1]"
                    )
                betas.append(v)
            if len(betas) != len(samples):
                raise ValueError(
                    f"{path}: probe {probe_id} has {len(betas)} values, "
                    f"expected {len(samples)} matching header columns"
                )
            probe_betas[probe_id] = betas

    return samples, probe_betas


# ---------- internals ----------


def _open_text(path: str | Path, mode: str) -> IO[str]:
    """Open a possibly-gzipped text file."""

    p = Path(path)
    if p.suffix == ".gz":
        return _io.TextIOWrapper(gzip.open(p, mode.replace("t", "b")), encoding="utf-8")
    return open(p, mode, encoding="utf-8")


def write_jsonl_atomic(path: str | Path, records: Iterable[BaseModel]) -> int:
    """Atomic JSONL writer: write to a sibling tempfile, then `os.replace`.

    The temp filename inserts `.tmp` BEFORE the final suffix
    (`catalog.jsonl.gz` → `catalog.jsonl.tmp.gz`) so `_open_text` still
    dispatches to gzip when the user requested a `.gz` output. Earlier
    naïve `with_suffix(suffix + ".tmp")` produced `catalog.jsonl.gz.tmp`,
    which `_open_text` saw as a plain-text suffix and wrote uncompressed —
    leaving an unreadable file with a `.gz` extension after the rename.
    """

    import os

    p = Path(path)
    tmp = p.with_name(f"{p.stem}.tmp{p.suffix}")
    n = write_jsonl(tmp, records)
    os.replace(tmp, p)
    return n


def write_json(path: str | Path, payload: object) -> None:
    with _open_text(path, "wt") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def read_sample_subtypes(path: str | Path) -> dict[str, str]:
    """Load a sample → subtype mapping from a TSV with columns `sample_id, subtype`.

    Used by the V2 subtype-aware backend factory to split a single beta matrix
    into per-subtype submatrices (e.g. PAM50 LumA / LumB / Basal / HER2-enriched).

    Conflicting duplicates are an error: a sample mapped to two different
    subtypes silently routed into the wrong cohort under the previous
    last-wins behavior. Identical duplicate rows (same sample, same subtype)
    are tolerated since they're harmless.
    """

    out: dict[str, str] = {}
    for row in read_tsv(path):
        try:
            sid, subtype = row["sample_id"], row["subtype"]
        except KeyError as e:
            raise ValueError(
                f"{path}: expected columns sample_id, subtype — missing {e}"
            ) from e
        if sid in out and out[sid] != subtype:
            raise ValueError(
                f"{path}: sample {sid!r} mapped to conflicting subtypes "
                f"{out[sid]!r} and {subtype!r}"
            )
        out[sid] = subtype
    return out
