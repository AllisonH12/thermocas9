# Roth follow-up — native EPIC v2 confirmation

**Status.** Draft. Edit before sending.
**Audience.** Dr. Luca Roth (and the corresponding-author team on Roth et al., *Nature* 2026).
**Framework.** Code/results at `memo-2026-04-22-z` of <https://github.com/AllisonH12/thermocas9>.
**Purpose.** Close the loop on the HM450-intersect shortcut from our earlier note, confirm the native EPIC v2 path reproduces the same biological story, and offer the annotated shortlist for their own sanity-checking. No ask.

---

## Subject line options

- `ThermoCas9 target-ranking: native EPIC v2 confirmation on GSE322563`
- `Follow-up on GSE322563: native EPIC v2 run matches the HM450-intersect result`

---

## Draft body

Dr. Roth,

Following up on our earlier note about GSE322563. A correction and a confirmation:

**Correction.** Our first message cited the accession as `GSE32256`; the correct accession is `GSE322563` as reported in your supplementary. Apologies for the typo; it delayed our ingest by about a day.

**Confirmation.** The first draft of our pipeline harmonized EPIC v2 probe IDs to HM450 by stripping the `_BC##`/`_TC##`/`_TO##`/`_BO##` beadchip-design suffix and intersecting with the HM450 universe (~80.7 % retention). We've since re-run the analysis against the **native EPIC v2 probe set** — 147,928 probes on chr5/6/10, lifted hg38 → hg19 from GPL33022 — with no HM450 intersect.

The two paths agree on the headline conclusions:

1. **Δβ distributions** on the MCF-7 vs MCF-10A differential are consistent between HM450-intersect and native runs.
2. **Validated-target AUC agrees to within the n=2/2 noise floor.** V2.5 AUC on the n = 3 Roth-validated probes is 0.990 (HM450-intersect) vs 0.986 (native EPIC v2) — a 0.003 movement, well inside the tied-band noise floor on this cohort. Per-positive ranks shift more (e.g. ESR1: 2,575 → 5,746; EGFLAM: 64,433 → 158,502; full table at `examples/validated_positive_ranks.md`), reflecting the larger native catalog (~5.2M vs ~3.0M candidates), but the AUC ranking conclusion is preserved.
3. **Top-K windows share several high-scoring examples** between the two paths (e.g. promoter CpGs at *KCNIP2*, *CALHM2*, *CELF2* surface in both top-20 windows), but the windows are *not* position-preserved — both top-20s are 20-record windows inside several-hundred-record tied bands at K=100, so the absolute ordering inside each window depends on the documented `candidate_id` ascending tie-break.
4. **The HM450-intersect shortcut did not distort the V2.5 scoring claim.** We had to check: given the retention rate, it was a reasonable reviewer question.

We take this as closing the obvious "did you lose something by intersecting to HM450?" concern. The native path is now the primary in our memo (§4.4, §5.2).

**Two things that may be useful to you, separate from the method claim**:

- **Annotation layer.** For each top hit we now record nearest gene, TSS distance, promoter / gene-body / intergenic classification, CpG-island context (island / shore / shelf / open-sea), RepeatMasker overlap, and ENCODE DNase-HS cluster breadth. On the MCF-7 shortlist, roughly 4 of 20 hits sit in repeat-rich regions (primarily L1 / ERVL-MaLR elements) or have no DNase-cluster support — flagged in the companion .md, not re-ranked away.
- **Cohort n = 2 caveat.** GSE322563 is two tumor + two normal replicates, so `p_trust` saturates low across the shortlist. We surface this as a per-candidate CAUTION flag ("sparse evidence") rather than hiding it. The V2.5 numbers are meaningful at the ranking level but should not be read as per-site significance.

If it's useful, we're happy to share the scored JSONLs and the annotated shortlist TSV + Markdown companion. Everything is reproducible from the tagged revision above; the build commands are in `examples/gse322563_native_roth_labels/`.

No action required. We just wanted you to know the EPIC v2 ingest is working and doesn't change the biological story you reported.

Best,
Allison Huang
Columbia University
<allisonhmercer@gmail.com>

---

## Notes for the sender

- **Tone**: friendly, technical, no method-evangelism. The Roth group shipped the paper that made this work possible; this note is a courtesy and a sanity-check offer, not a sell.
- **Don't oversell V2.5**: the message leads with the confirmation that the HM450 shortcut wasn't hiding anything, not with "we built a better scorer." V2.5 gets one parenthetical mention.
- **Don't attach the repo tarball**. Link to the tagged GitHub revision. If they want to run anything they can clone.
- **Include both the TSV and the .md** if they ask for the shortlist. The .md is the one they'll actually read.
- **Revision checkpoint**: if this note sits more than 48 h before sending, re-verify the tag reference (`git rev-parse memo-2026-04-22-z` should still resolve).
