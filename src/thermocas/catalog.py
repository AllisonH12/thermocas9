"""Genome-wide PAM catalog builder.

Streams a FASTA reference chromosome by chromosome, runs the configured PAM
model over each sequence, and emits `CandidateSite` records. Memory cost
scales with the longest chromosome, not the genome.

The output is suitable for direct JSONL serialization via `io.write_jsonl`.
The catalog is cohort-independent and gets reused across cohorts.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from pathlib import Path
from typing import NamedTuple

from bisect import bisect_left

from thermocas.io import iter_fasta
from thermocas.models import CandidateSite, Strand
from thermocas.pam_model import PamModel, find_pam_matches, reverse_complement

#: Width of the local sequence context attached to each CandidateSite (±50 bp).
LOCAL_CONTEXT_HALF_WIDTH = 50


def _local_context(seq: str, n: int, critical_c_pos: int, strand: Strand) -> str:
    """Extract the ±50 bp window around `critical_c_pos`, oriented on `strand`.

    The `CandidateSite` contract states that `local_seq_100bp` is on the
    *indicated strand*. For minus-strand candidates we reverse-complement the
    forward window so the PAM in `match.sequence` actually appears in the
    returned context, oriented 5'→3' on the minus strand.
    """

    lo = max(0, critical_c_pos - LOCAL_CONTEXT_HALF_WIDTH)
    hi = min(n, critical_c_pos + LOCAL_CONTEXT_HALF_WIDTH)
    fwd_window = seq[lo:hi]
    if strand == Strand.MINUS:
        return reverse_complement(fwd_window)
    return fwd_window

# A region filter: returns True if a candidate at (chrom, critical_c_pos) should
# be kept. Use to restrict to e.g. a gene-body BED.
RegionFilter = Callable[[str, int], bool]


def probe_window_filter(
    probe_annotation_path: str | Path, window_bp: int
) -> RegionFilter:
    """Build a `RegionFilter` that keeps only candidates within ±`window_bp`
    of any probe in `probe_annotation_path` (a `probe_id, chrom, pos` TSV).

    Reduces a genome-wide ThermoCas9 catalog from ~370M candidates to ~10M
    by dropping everything that would later score as `UNOBSERVED` anyway.
    Catalog build time and on-disk size both drop ~95%.
    """

    from thermocas.io import read_tsv  # local import keeps catalog.py import-light

    sorted_positions: dict[str, list[int]] = {}
    for row in read_tsv(probe_annotation_path):
        try:
            chrom = row["chrom"]
            pos = int(row["pos"])
        except (KeyError, ValueError) as e:
            raise ValueError(
                f"{probe_annotation_path}: probe rows must have chrom + pos: {e}"
            ) from e
        sorted_positions.setdefault(chrom, []).append(pos)
    for chrom in sorted_positions:
        sorted_positions[chrom].sort()

    def keep(chrom: str, critical_c_pos: int) -> bool:
        positions = sorted_positions.get(chrom)
        if not positions:
            return False
        idx = bisect_left(positions, critical_c_pos)
        candidates_to_check = []
        if idx < len(positions):
            candidates_to_check.append(positions[idx])
        if idx > 0:
            candidates_to_check.append(positions[idx - 1])
        return any(abs(p - critical_c_pos) <= window_bp for p in candidates_to_check)

    return keep


class CatalogStats(NamedTuple):
    """Summary returned by `build_catalog` for logging / sanity-checking."""

    n_chromosomes: int
    n_candidates: int
    n_filtered_out: int


def build_catalog(
    fasta_path: str | Path,
    pam_model: PamModel,
    region_filter: RegionFilter | None = None,
    annotator: Callable[[CandidateSite], CandidateSite] | None = None,
) -> tuple[list[CandidateSite], CatalogStats]:
    """Scan a reference FASTA and return all matching CandidateSites.

    Args:
        fasta_path: path to a FASTA (optionally `.gz`).
        pam_model: configured `PamModel`.
        region_filter: optional `(chrom, critical_c_pos) -> bool`. False drops the candidate.
        annotator: optional final-pass annotator, e.g. attach `nearest_gene`.

    Returns:
        (candidates, stats). Candidates carry deterministic `candidate_id`s of
        the form `chrom:critical_c_pos:strand:family[:idx]`.
    """

    candidates: list[CandidateSite] = []
    n_chromosomes = 0
    n_filtered_out = 0
    seen_ids: set[str] = set()

    for chrom, seq in iter_fasta(fasta_path):
        n_chromosomes += 1
        seq = seq.upper()
        n = len(seq)
        for match in find_pam_matches(seq, pam_model):
            if region_filter is not None and not region_filter(chrom, match.critical_c_pos):
                n_filtered_out += 1
                continue

            base_id = f"{chrom}:{match.critical_c_pos}{match.strand.value}:{match.family}"
            cid = base_id
            # extremely unlikely with current PAM model + strand keys, but harden anyway
            suffix = 1
            while cid in seen_ids:
                cid = f"{base_id}:{suffix}"
                suffix += 1
            seen_ids.add(cid)

            local = _local_context(seq, n, match.critical_c_pos, match.strand)

            site = CandidateSite(
                candidate_id=cid,
                chrom=chrom,
                critical_c_pos=match.critical_c_pos,
                strand=match.strand,
                pam=match.sequence,
                pam_family=match.family,
                is_cpg_pam=match.is_cpg,
                local_seq_100bp=local,
            )
            if annotator is not None:
                site = annotator(site)
            candidates.append(site)

    return candidates, CatalogStats(
        n_chromosomes=n_chromosomes,
        n_candidates=len(candidates),
        n_filtered_out=n_filtered_out,
    )


def stream_catalog(
    fasta_path: str | Path,
    pam_model: PamModel,
    region_filter: RegionFilter | None = None,
) -> Iterator[CandidateSite]:
    """Same as `build_catalog` but yields candidates lazily, one chromosome at a time.

    Use for genome-scale builds where holding the entire candidate list in memory
    is undesirable. Pair with `io.write_jsonl` to stream straight to disk.
    """

    seen_ids: set[str] = set()
    for chrom, seq in iter_fasta(fasta_path):
        seq = seq.upper()
        n = len(seq)
        for match in find_pam_matches(seq, pam_model):
            if region_filter is not None and not region_filter(chrom, match.critical_c_pos):
                continue
            base_id = f"{chrom}:{match.critical_c_pos}{match.strand.value}:{match.family}"
            cid = base_id
            suffix = 1
            while cid in seen_ids:
                cid = f"{base_id}:{suffix}"
                suffix += 1
            seen_ids.add(cid)

            yield CandidateSite(
                candidate_id=cid,
                chrom=chrom,
                critical_c_pos=match.critical_c_pos,
                strand=match.strand,
                pam=match.sequence,
                pam_family=match.family,
                is_cpg_pam=match.is_cpg,
                local_seq_100bp=_local_context(seq, n, match.critical_c_pos, match.strand),
            )
