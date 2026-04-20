"""PAM model loader and matcher.

Reads a YAML PAM model and returns concrete `PamMatch` objects for every
occurrence of a configured PAM family in a forward sequence string. Reverse
strand handled by reverse-complementing and re-scanning, so the caller gets
strand-aware matches in one call to `find_pam_matches`.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

import yaml
from pydantic import BaseModel

from thermocas.models import PamFamily, PamMatch, Strand

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a nucleotide sequence (preserves N)."""

    return seq.translate(_COMPLEMENT)[::-1]


class PamModel(BaseModel):
    """Container for a configured set of ThermoCas9 PAM families."""

    pam_families: list[PamFamily]

    @classmethod
    def from_yaml(cls, path: str | Path) -> PamModel:
        with Path(path).open() as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)

    def get(self, name: str) -> PamFamily:
        for fam in self.pam_families:
            if fam.name == name:
                return fam
        raise KeyError(f"Unknown PAM family: {name}")


def _scan_forward(seq: str, family: PamFamily, strand: Strand) -> Iterator[PamMatch]:
    """Yield PamMatch objects for a single strand, including overlapping matches.

    `re.finditer` advances past each match's full span and so silently drops
    overlapping motifs (`AAAACCACCA` matches `NNNNCCA` at both offset 0 and 3,
    but a naive scan returns only offset 0). For genome-wide PAM enumeration
    this is a correctness bug — every valid window must be reported. The
    lookahead `(?=(...))` lets `finditer` advance one base at a time while
    still capturing the full PAM via group 1.
    """

    pattern = re.compile(f"(?=({family.regex}))")
    for m in pattern.finditer(seq):
        match_seq = m.group(1)
        start = m.start()
        critical_c = start + family.critical_c_offset
        yield PamMatch(
            family=family.name,
            sequence=match_seq,
            start=start,
            strand=strand,
            critical_c_pos=critical_c,
            is_cpg=family.is_cpg,
            weight=family.weight,
        )


def find_pam_matches(seq: str, model: PamModel) -> list[PamMatch]:
    """Return every PAM occurrence on either strand of `seq`.

    Coordinates are reported in the forward-strand reference frame even for
    minus-strand hits — i.e. `critical_c_pos` is always the index in the input
    `seq`, not in the reverse-complement string.
    """

    seq = seq.upper()
    n = len(seq)
    matches: list[PamMatch] = []
    for fam in model.pam_families:
        # forward strand
        matches.extend(_scan_forward(seq, fam, Strand.PLUS))
        # reverse strand: scan revcomp, then translate coords back
        rc = reverse_complement(seq)
        for m in _scan_forward(rc, fam, Strand.MINUS):
            # critical C in revcomp coords → forward coord
            #   forward_pos = (n - 1) - revcomp_pos
            fwd_critical = (n - 1) - m.critical_c_pos
            fwd_start = (n - 1) - (m.start + len(m.sequence) - 1)
            matches.append(
                PamMatch(
                    family=fam.name,
                    sequence=m.sequence,  # PAM as written on the minus strand (revcomp space)
                    start=fwd_start,
                    strand=Strand.MINUS,
                    critical_c_pos=fwd_critical,
                    is_cpg=fam.is_cpg,
                    weight=fam.weight,
                )
            )
    matches.sort(key=lambda m: (m.critical_c_pos, m.strand.value, m.family))
    return matches
