"""Probe-to-site evidence classification.

Methylation arrays measure a predefined subset of CpG sites — HM450 covers
~450,000 probe-associated CpGs. The framework treats this sparse coverage as
a first-class problem: every candidate carries an `EvidenceClass` indicating
how directly the assay observation maps onto the candidate's critical PAM
cytosine.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from thermocas.models import EvidenceClass, EvidenceThresholds


@dataclass(frozen=True)
class ProbeRecord:
    """A minimal probe descriptor — what every methylation backend must produce.

    Backend implementations (GDC, raw IDAT, bisulfite) produce these.
    Coordinates are 0-indexed forward-strand positions of the CpG cytosine.
    """

    probe_id: str
    chrom: str
    pos: int


def classify_evidence(
    distance_bp: int | None, thresholds: EvidenceThresholds
) -> EvidenceClass:
    """Map an absolute distance to the nearest informative probe → EvidenceClass.

    `None` distance means no probe was found within `regional_bp` → UNOBSERVED.
    """

    if distance_bp is None or distance_bp > thresholds.regional_bp:
        return EvidenceClass.UNOBSERVED
    if distance_bp <= thresholds.exact_bp:
        return EvidenceClass.EXACT
    if distance_bp <= thresholds.proximal_close_bp:
        return EvidenceClass.PROXIMAL_CLOSE
    if distance_bp <= thresholds.proximal_bp:
        return EvidenceClass.PROXIMAL
    return EvidenceClass.REGIONAL


class EvidenceClassifier:
    """Find the nearest probe for each candidate critical-C position.

    Indexes probes per chromosome in sorted order so each lookup is O(log n).
    """

    def __init__(self, probes: Iterable[ProbeRecord], thresholds: EvidenceThresholds):
        self._thresholds = thresholds
        # chrom → sorted list of (pos, probe)
        self._by_chrom: dict[str, list[tuple[int, ProbeRecord]]] = {}
        for p in probes:
            self._by_chrom.setdefault(p.chrom, []).append((p.pos, p))
        for lst in self._by_chrom.values():
            lst.sort(key=lambda t: t[0])

    def nearest(self, chrom: str, pos: int) -> tuple[ProbeRecord | None, int | None]:
        """Return (nearest probe, |distance|) on `chrom` to `pos`. None when no probes."""

        rows = self._by_chrom.get(chrom)
        if not rows:
            return None, None

        # binary search
        lo, hi = 0, len(rows) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if rows[mid][0] < pos:
                lo = mid + 1
            else:
                hi = mid
        # check `lo` and `lo - 1`
        candidates: list[tuple[int, ProbeRecord]] = []
        if 0 <= lo < len(rows):
            candidates.append(rows[lo])
        if lo - 1 >= 0:
            candidates.append(rows[lo - 1])
        best_pos, best_probe = min(candidates, key=lambda t: abs(t[0] - pos))
        return best_probe, abs(best_pos - pos)

    def classify(self, chrom: str, pos: int) -> tuple[EvidenceClass, ProbeRecord | None, int | None]:
        """Return (evidence class, nearest probe, distance) for one critical-C position."""

        probe, dist = self.nearest(chrom, pos)
        ec = classify_evidence(dist, self._thresholds)
        if ec == EvidenceClass.UNOBSERVED:
            return ec, None, None
        return ec, probe, dist
