"""V3 — gRNA protospacer scoring.

The framework's V1/V2 scoring ranks PAM sites by methylation-driven cancer
selectivity. None of that says anything about whether the resulting 20-nt
protospacer (the actual gRNA) is well-designed. Two sites can have identical
methylation evidence but very different sgRNA quality (one has a 7-base
poly-T run, the other doesn't). V3 closes that gap.

Scoring is intentionally heuristic — a few well-known sgRNA design rules
combined into a single `SpacerScore`. The scores are interpretable individually
so callers can audit which dimension dominated. A future PR can swap in a
learned model (e.g. Doench Rule Set 2) by replacing `score_spacer`.

Components (each in [0, 1], higher is better):
    * gc_content_score : favors 40–60% GC, penalizes <30% or >70%
    * tm_score         : favors Tm in 60–75 °C
    * runs_score       : penalizes mononucleotide runs ≥ 5
    * hairpin_score    : penalizes 4-bp self-complementary stretches
                         (rough RNA hairpin proxy)

`final_score` is the geometric mean of the four components — one bad component
pulls the spacer down rather than averaging away.
"""

from __future__ import annotations

from thermocas.models import CandidateSite, SpacerScore, Strand
from thermocas.pam_model import reverse_complement

#: Length of the protospacer (Type II Cas9 = 20 nt).
SPACER_LEN = 20


# ---------- public API ----------


def extract_spacer(candidate: CandidateSite) -> str | None:
    """Pull the 20-nt protospacer immediately upstream of the candidate's PAM.

    `local_seq_100bp` is stored on the candidate's strand (V1 strand-orientation
    fix), with the critical PAM cytosine roughly in the middle. The PAM sits at
    `[half - critical_c_offset, half + (pam_width - critical_c_offset)]`; the
    protospacer is the 20 nt **before** that, on the same strand.

    Returns None when there isn't enough flanking context (chrom edges).
    """

    seq = candidate.local_seq_100bp
    if not seq:
        return None
    # Locate the PAM in the local context. PAM is on the indicated strand, so
    # local_seq is already oriented; search for the literal PAM string.
    pam_start = seq.find(candidate.pam)
    if pam_start < 0 or pam_start < SPACER_LEN:
        return None
    spacer = seq[pam_start - SPACER_LEN : pam_start]
    if len(spacer) != SPACER_LEN:
        return None
    if "N" in spacer:
        return None
    return spacer


def score_spacer(
    candidate: CandidateSite,
    spacer: str | None = None,
) -> SpacerScore | None:
    """Compute a `SpacerScore` for a candidate. Returns None when no spacer
    can be extracted (chrom edge, or `local_seq_100bp` doesn't include the PAM)."""

    if spacer is None:
        spacer = extract_spacer(candidate)
    if spacer is None or len(spacer) != SPACER_LEN:
        return None

    gc_frac = _gc_fraction(spacer)
    tm = _melting_temp(spacer)
    longest_run = _longest_mononucleotide_run(spacer)

    gc_score = _gc_content_score(gc_frac)
    tm_score = _tm_score(tm)
    runs_score = _runs_score(longest_run)
    hairpin_score = _hairpin_score(spacer)

    final = _geometric_mean(gc_score, tm_score, runs_score, hairpin_score)

    return SpacerScore(
        candidate_id=candidate.candidate_id,
        spacer_seq=spacer,
        gc_fraction=gc_frac,
        tm_C=tm,
        longest_run=longest_run,
        gc_content_score=gc_score,
        tm_score=tm_score,
        runs_score=runs_score,
        hairpin_score=hairpin_score,
        final_score=final,
    )


# ---------- component scorers ----------


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in "GC")
    return gc / len(seq)


def _melting_temp(seq: str) -> float:
    """Wallace-rule Tm in °C: 2·(A+T) + 4·(G+C). Cheap and adequate for the
    short, fixed-length protospacer; not meant to compete with nearest-neighbor
    thermodynamic estimates."""

    a_t = sum(1 for c in seq if c in "AT")
    g_c = sum(1 for c in seq if c in "GC")
    return 2.0 * a_t + 4.0 * g_c


def _longest_mononucleotide_run(seq: str) -> int:
    if not seq:
        return 0
    longest = 1
    current = 1
    for prev, this in zip(seq, seq[1:], strict=False):
        if prev == this:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def _gc_content_score(gc_frac: float) -> float:
    """Triangular function: peak (1.0) in [0.40, 0.60], 0.0 outside [0.20, 0.80]."""
    if 0.40 <= gc_frac <= 0.60:
        return 1.0
    if gc_frac < 0.20 or gc_frac > 0.80:
        return 0.0
    if gc_frac < 0.40:
        return (gc_frac - 0.20) / 0.20
    return (0.80 - gc_frac) / 0.20  # gc_frac > 0.60


def _tm_score(tm_c: float) -> float:
    """Triangular: peak in [60, 75] °C, 0 outside [45, 90]."""
    if 60.0 <= tm_c <= 75.0:
        return 1.0
    if tm_c < 45.0 or tm_c > 90.0:
        return 0.0
    if tm_c < 60.0:
        return (tm_c - 45.0) / 15.0
    return (90.0 - tm_c) / 15.0


def _runs_score(longest_run: int) -> float:
    """1.0 for runs ≤ 4; linearly down to 0 at runs ≥ 8."""
    if longest_run <= 4:
        return 1.0
    if longest_run >= 8:
        return 0.0
    return (8 - longest_run) / 4.0


def _hairpin_score(seq: str) -> float:
    """Crude self-complementarity penalty.

    Slides a 4-mer window across the spacer and checks whether its
    reverse-complement appears later in the same spacer (≥ 6 bp downstream).
    Each unique hit is a penalty; capped at 4 hits so single bad spacers
    saturate to 0.0 rather than going negative.
    """

    n = len(seq)
    if n < 8:
        return 1.0
    hits = 0
    seen: set[str] = set()
    for i in range(n - 4):
        kmer = seq[i : i + 4]
        if "N" in kmer or kmer in seen:
            continue
        seen.add(kmer)
        rc = reverse_complement(kmer)
        # require the complement window to start ≥ 6 bp after `i` so we don't
        # count immediate adjacencies that aren't real hairpins
        if rc in seq[i + 6 :]:
            hits += 1
    capped = min(hits, 4)
    return 1.0 - capped / 4.0


def _geometric_mean(*xs: float) -> float:
    if not xs:
        return 0.0
    if any(x <= 0.0 for x in xs):
        return 0.0
    product = 1.0
    for x in xs:
        product *= x
    return product ** (1.0 / len(xs))


__all__ = [
    "SPACER_LEN",
    "extract_spacer",
    "score_spacer",
]
