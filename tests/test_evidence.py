"""Evidence classifier tests."""

from __future__ import annotations

from thermocas.evidence import EvidenceClassifier, ProbeRecord, classify_evidence
from thermocas.models import EvidenceClass, EvidenceThresholds


def test_classify_evidence_bins():
    t = EvidenceThresholds()
    assert classify_evidence(0, t) == EvidenceClass.EXACT
    assert classify_evidence(t.exact_bp, t) == EvidenceClass.EXACT
    assert classify_evidence(t.exact_bp + 1, t) == EvidenceClass.PROXIMAL_CLOSE
    assert classify_evidence(t.proximal_close_bp, t) == EvidenceClass.PROXIMAL_CLOSE
    assert classify_evidence(t.proximal_close_bp + 1, t) == EvidenceClass.PROXIMAL
    assert classify_evidence(t.proximal_bp, t) == EvidenceClass.PROXIMAL
    assert classify_evidence(t.proximal_bp + 1, t) == EvidenceClass.REGIONAL
    assert classify_evidence(t.regional_bp, t) == EvidenceClass.REGIONAL
    assert classify_evidence(t.regional_bp + 1, t) == EvidenceClass.UNOBSERVED
    assert classify_evidence(None, t) == EvidenceClass.UNOBSERVED


def test_classifier_finds_nearest_probe():
    probes = [
        ProbeRecord(probe_id="cg001", chrom="chr1", pos=100),
        ProbeRecord(probe_id="cg002", chrom="chr1", pos=200),
        ProbeRecord(probe_id="cg003", chrom="chr1", pos=10_000),
        ProbeRecord(probe_id="cg004", chrom="chr2", pos=500),
    ]
    cls = EvidenceClassifier(probes, EvidenceThresholds())

    # exact hit
    ec, probe, dist = cls.classify("chr1", 100)
    assert ec == EvidenceClass.EXACT
    assert probe is not None and probe.probe_id == "cg001"
    assert dist == 0

    # proximal_close
    ec, probe, dist = cls.classify("chr1", 220)
    assert ec == EvidenceClass.PROXIMAL_CLOSE
    assert probe is not None and probe.probe_id == "cg002"
    assert dist == 20

    # proximal
    ec, probe, dist = cls.classify("chr1", 145)
    assert ec == EvidenceClass.PROXIMAL
    assert dist == 45

    # regional
    ec, probe, dist = cls.classify("chr1", 700)
    assert ec == EvidenceClass.REGIONAL
    assert dist == 500

    # unobserved
    ec, probe, dist = cls.classify("chr1", 50_000)
    assert ec == EvidenceClass.UNOBSERVED
    assert probe is None
    assert dist is None


def test_classifier_unknown_chrom():
    cls = EvidenceClassifier(
        [ProbeRecord(probe_id="cg001", chrom="chr1", pos=100)], EvidenceThresholds()
    )
    ec, probe, dist = cls.classify("chrX", 100)
    assert ec == EvidenceClass.UNOBSERVED
    assert probe is None
    assert dist is None
