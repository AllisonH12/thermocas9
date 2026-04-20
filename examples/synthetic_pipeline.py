"""End-to-end synthetic demo of the framework — runs without any real data.

Builds a small fake reference, scans for ThermoCas9 PAMs, fabricates probe
methylation values for a 'tumor' and a 'normal' cohort with a deliberate
selectivity differential, classifies evidence per candidate, and prints a
ranked candidate table.

Run:
    python examples/synthetic_pipeline.py
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

# allow `python examples/synthetic_pipeline.py` from repo root without install
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from thermocas.config import load_cohort_config  # noqa: E402
from thermocas.evidence import EvidenceClassifier, ProbeRecord  # noqa: E402
from thermocas.models import (  # noqa: E402
    CandidateSite,
    EvidenceClass,
    MethylationObservation,
)
from thermocas.pam_model import PamModel, find_pam_matches, reverse_complement  # noqa: E402
from thermocas.scoring import score_candidate  # noqa: E402

CONFIG_DIR = REPO_ROOT / "config"

# ---------- 1. Fake reference ----------

rng = random.Random(42)
SEQ_LEN = 5_000
SEQ = "".join(rng.choices("ACGT", k=SEQ_LEN))
CHROM = "synth1"

# ---------- 2. PAM scan ----------

pam_model = PamModel.from_yaml(CONFIG_DIR / "pam_model.yaml")
matches = find_pam_matches(SEQ, pam_model)
print(f"Found {len(matches)} PAM hits in {SEQ_LEN}-bp synthetic chromosome.")

candidates: list[CandidateSite] = []
for i, m in enumerate(matches):
    fwd_window = SEQ[max(0, m.critical_c_pos - 50): m.critical_c_pos + 50]
    # local_seq is on the *indicated strand* — revcomp for minus-strand hits so
    # the PAM in `m.sequence` actually appears in the returned context.
    local = reverse_complement(fwd_window) if m.strand.value == "-" else fwd_window
    candidates.append(
        CandidateSite(
            candidate_id=f"{CHROM}:{m.critical_c_pos}{m.strand.value}:{m.family}:{i}",
            chrom=CHROM,
            critical_c_pos=m.critical_c_pos,
            strand=m.strand,
            pam=m.sequence,
            pam_family=m.family,
            is_cpg_pam=m.is_cpg,
            local_seq_100bp=local,
            nearest_gene="GENE_FAKE",
            regulatory_context="promoter" if m.critical_c_pos < SEQ_LEN // 2 else "gene_body",
        )
    )

# ---------- 3. Fake methylation probes ----------

# Place probes randomly; some will hit candidates exactly, some won't.
probes = [
    ProbeRecord(probe_id=f"cg{p:06d}", chrom=CHROM, pos=p)
    for p in sorted(rng.sample(range(SEQ_LEN), k=200))
]

# ---------- 4. Cohort config + evidence classification ----------

cohort = load_cohort_config(CONFIG_DIR / "cohorts" / "brca_example.yaml")
classifier = EvidenceClassifier(probes, cohort.evidence_thresholds)


def fake_betas(critical_c_pos: int) -> tuple[float, float]:
    """Half the candidates are 'tumor-hypomethylated' (selective), half are not."""

    if critical_c_pos % 2 == 0:
        return 0.05, 0.85   # selective: low tumor, high normal
    return 0.55, 0.55       # non-selective: equal


# ---------- 5. Build observations ----------

observations: list[MethylationObservation] = []
for cand in candidates:
    ec, probe, dist = classifier.classify(cand.chrom, cand.critical_c_pos)
    if ec == EvidenceClass.UNOBSERVED:
        observations.append(
            MethylationObservation(
                candidate_id=cand.candidate_id,
                cohort_name=cohort.name,
                evidence_class=ec,
            )
        )
        continue
    bt, bn = fake_betas(cand.critical_c_pos)
    observations.append(
        MethylationObservation(
            candidate_id=cand.candidate_id,
            cohort_name=cohort.name,
            evidence_class=ec,
            evidence_distance_bp=dist,
            probe_id=probe.probe_id if probe else None,
            beta_tumor_mean=bt,
            beta_tumor_q25=max(0.0, bt - 0.05),
            beta_tumor_q75=min(1.0, bt + 0.05),
            n_samples_tumor=400,
            beta_normal_mean=bn,
            beta_normal_q25=max(0.0, bn - 0.07),
            beta_normal_q75=min(1.0, bn + 0.07),
            n_samples_normal=80,
        )
    )

# ---------- 6. Score and rank ----------

scored = []
for cand, obs in zip(candidates, observations, strict=True):
    fam = pam_model.get(cand.pam_family)
    scored.append(score_candidate(cand, obs, fam, cohort))

scored.sort(key=lambda s: s.final_score, reverse=True)

# ---------- 7. Print top + bottom ----------

print()
print(f"Top 5 ranked candidates (cohort: {cohort.name}):")
print("-" * 100)
print(f"{'rank':>4}  {'candidate_id':<48}  {'evidence':<14}  "
      f"{'tumor_beta':>10}  {'normal_beta':>11}  {'score':>8}")
for i, s in enumerate(scored[:5], 1):
    bt = s.observation.beta_tumor_mean
    bn = s.observation.beta_normal_mean
    print(
        f"{i:>4}  {s.candidate.candidate_id:<48}  "
        f"{s.observation.evidence_class.value:<14}  "
        f"{(f'{bt:.2f}' if bt is not None else '   --'):>10}  "
        f"{(f'{bn:.2f}' if bn is not None else '    --'):>11}  "
        f"{s.final_score:>8.3f}"
    )

print()
print(f"Bottom 5 ranked candidates (cohort: {cohort.name}):")
print("-" * 100)
for i, s in enumerate(scored[-5:], len(scored) - 4):
    bt = s.observation.beta_tumor_mean
    bn = s.observation.beta_normal_mean
    print(
        f"{i:>4}  {s.candidate.candidate_id:<48}  "
        f"{s.observation.evidence_class.value:<14}  "
        f"{(f'{bt:.2f}' if bt is not None else '   --'):>10}  "
        f"{(f'{bn:.2f}' if bn is not None else '    --'):>11}  "
        f"{s.final_score:>8.3f}"
    )

# ---------- 8. Summary stats ----------

n_total = len(scored)
n_observed = sum(1 for s in scored if s.observation.evidence_class != EvidenceClass.UNOBSERVED)
n_positive = sum(1 for s in scored if s.final_score > 0)
print()
print(f"Summary: {n_total} candidates · {n_observed} observed ({n_observed/n_total:.0%}) · "
      f"{n_positive} with final_score > 0 ({n_positive/n_total:.0%})")
