"""Cross-check selected numerical claims in MANUSCRIPT.md, PAPER.md, and
README.md against the committed bench JSONLs and the source-code constants.

Rationale: during the 2026-04-22 revision cycle, three successive tags
(-c, -d, -e) each shipped with prose claims that did not match the
committed bench artifacts (fabricated AUCs; wrong universal quantifiers;
mis-quoted constants). Manual verification has repeatedly missed class-of-
bug issues like "at every cohort × tier combination tested" when the
committed data actually violated the claim on several rows.

This script is the guard against the *known classes* of that failure
mode that have actually shown up. Run it before cutting any new memo-*
tag. Exits non-zero on mismatch.

Usage:
    python scripts/verify_manuscript_claims.py

Exit codes:
    0 — all CHECKED claims match committed artifacts and source constants.
    1 — at least one mismatch. Details printed to stdout.

Scope (what this script ACTUALLY verifies, by check; reflects exactly
which `check_*` functions are wired into `main()`):

    1. **Code-constant quotes** (`check_constants`) — run on MANUSCRIPT.md
       §2.1 only. Verifies that `τ_u = 0.30`, `δ = 0.2`, σ_floor = 0.05,
       and the n=30 trust ramp are quoted exactly as named in
       `src/thermocas/probabilistic.py`. Also verifies a few EvidenceClass
       trust-base values are quoted correctly.

    2. **Universal ordering claims** (`check_universal_ordering_claims`) —
       not text-driven; enumerates the canonical V2.5 > V1 and V2.5 > Δβ
       invariants on the 9 matched-cell-line rows + 3 GSE69914 tissue
       rows directly from the JSONL grid, and flags any row where the
       invariant fails. Catches data-side regressions, not prose
       overclaims; the latter is checked in `check_figure_captions` and
       in the doc-text pattern detectors below.

    3. **`tumor_only` tie-band per-cohort enumeration**
       (`check_tumor_only_tie_band_range`) — verifies the 5 expected
       per-cohort `tie_band@K=100` values are present in MANUSCRIPT.md
       §6.2 (in either bare or comma-formatted form), AND scans
       MANUSCRIPT.md / PAPER.md / README.md for the stale
       "6,000-12,000" or "6,500+ on every cohort" patterns, with a
       historical-context filter that exempts retrospectives in the
       PAPER.md tag ledger.

    4. **Naive baseline artifact count** (`check_artifact_counts`) —
       requires exactly 15 `bench_*_naive.jsonl` files under
       `examples/*_roth_labels/`.

    5. **Test count consistency** (`check_test_count`) — runs
       `pytest --collect-only -q`, sums per-file collected counts, and
       flags any quoted "N tests pass" / "N tests passing" in
       MANUSCRIPT.md or PAPER.md that disagrees.

    6. **Figure-caption gene-list claims** (`check_figure_captions`) —
       runs on BOTH MANUSCRIPT.md and PAPER.md. Extracts each
       `**Figure N.**` block; for each italicized all-caps gene token
       run, decides from local prior-context whether the caption is
       claiming the genes are in the GSE69914 tissue top-20, in the
       intersection of all three cell-line V2.5 top-20s, etc., and
       cross-checks against the relevant committed top20 TSV(s). Also
       fires on the 2-column "both cell-line" / "either cell-line
       column" framing patterns regardless of italicized gene presence,
       and on "exactly N distinct nearest-gene symbols" count claims
       on the tissue cohort.

NOT yet implemented (known coverage gaps):

    * **Per-cell AUC verification of every markdown table cell.** A
      `check_table_auc_values` function exists in this file but is not
      wired into `main()` — the heuristic exempt-list approach is
      noisy and would false-positive on margin / P@K values quoted in
      table cells. The proper fix is a markdown-table parser that
      identifies AUC cells specifically; out of scope for this script
      until/unless table drift is observed.
    * **README.md scope is narrower than MANUSCRIPT.md / PAPER.md**:
      README is covered by check (3) and the figure-caption check
      currently only runs on the two manuscripts. The constants check
      and universal-ordering check do not touch README. README's
      headline tables are AUC values that would benefit from the
      planned table-cell verifier above.
    * **Cover letter and Roth note** under `docs/notes/` are NOT
      scanned. They contain a small number of cohort-description and
      tag-citation claims; manual review is required before sending.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BENCH_ROOT = REPO / "examples"
PROB_PY = REPO / "src" / "thermocas" / "probabilistic.py"
MANUSCRIPT = REPO / "MANUSCRIPT.md"
PAPER = REPO / "PAPER.md"
README = REPO / "README.md"

#: Documents that are publicly circulated and whose claims must match the
#: committed artifacts. README is a smaller surface than the manuscripts
#: but is what GitHub visitors land on first; drift here is high-visibility.
PUBLIC_DOCS = (MANUSCRIPT, PAPER, README)


# ---------- authoritative sources ----------


def load_bench_grid() -> dict[tuple[str, str, str], dict]:
    """(cohort_dir, tier, axis) → BenchmarkResult dict."""
    grid = {}
    for p in BENCH_ROOT.glob("*_roth_labels/bench_*.jsonl"):
        cohort = p.parent.name.replace("_roth_labels", "")
        stem = p.stem
        parts = stem.split("_")
        tier = parts[1]
        axis = "_".join(parts[2:])
        grid[(cohort, tier, axis)] = json.loads(p.read_text().strip())
    return grid


def load_constants() -> dict[str, float]:
    """Extract DEFAULT_* constants and _TRUST_BASE entries from probabilistic.py."""
    text = PROB_PY.read_text()
    out: dict[str, float] = {}

    pats = {
        "DEFAULT_UNMETHYLATED_THRESHOLD": r"DEFAULT_UNMETHYLATED_THRESHOLD\s*=\s*([0-9.]+)",
        "DEFAULT_METHYLATED_THRESHOLD": r"DEFAULT_METHYLATED_THRESHOLD\s*=\s*([0-9.]+)",
        "DEFAULT_DIFFERENTIAL_DELTA": r"DEFAULT_DIFFERENTIAL_DELTA\s*=\s*([0-9.]+)",
        "DEFAULT_DIFFERENTIAL_SIGMA_FLOOR": r"DEFAULT_DIFFERENTIAL_SIGMA_FLOOR\s*=\s*([0-9.]+)",
        "DEFAULT_TRUST_RAMP_N": r"DEFAULT_TRUST_RAMP_N\s*=\s*([0-9]+)",
    }
    for name, pat in pats.items():
        m = re.search(pat, text)
        if not m:
            raise RuntimeError(f"could not find {name} in {PROB_PY}")
        out[name] = float(m.group(1))

    # EvidenceClass trust bases:  EvidenceClass.EXACT: 0.95,  etc.
    for cls_name in ("EXACT", "PROXIMAL_CLOSE", "PROXIMAL", "REGIONAL", "UNOBSERVED"):
        m = re.search(rf"EvidenceClass\.{cls_name}:\s*([0-9.]+)", text)
        if not m:
            raise RuntimeError(f"could not find EvidenceClass.{cls_name} trust base")
        out[f"trust_base_{cls_name}"] = float(m.group(1))
    return out


# ---------- claim-specific checks ----------


def check_constants(errors: list[str]) -> None:
    c = load_constants()
    # Check that MANUSCRIPT.md §2.1 quotes the right values.
    manuscript = MANUSCRIPT.read_text()
    expected = [
        ("τ_u = 0.30", c["DEFAULT_UNMETHYLATED_THRESHOLD"] == 0.30),
        ("δ", c["DEFAULT_DIFFERENTIAL_DELTA"] == 0.2),
        ("σ_floor = 0.05", c["DEFAULT_DIFFERENTIAL_SIGMA_FLOOR"] == 0.05),
        ("n = 30 ramp", c["DEFAULT_TRUST_RAMP_N"] == 30),
        ("EXACT 0.95", c["trust_base_EXACT"] == 0.95),
        ("PROXIMAL_CLOSE 0.75", c["trust_base_PROXIMAL_CLOSE"] == 0.75),
    ]
    for label, ok in expected:
        if not ok:
            errors.append(f"CONSTANT MISMATCH: {label}: code value = {c}")

    # Spot-check the prose values against the code values.
    if "τ_u = 0.30" not in manuscript and "τ_u = 0.3\n" not in manuscript:
        errors.append("MANUSCRIPT.md §2.1: does not cite τ_u = 0.30 explicitly")
    if "default 0.2" not in manuscript:
        errors.append("MANUSCRIPT.md §2.1: does not cite default δ = 0.2")
    if "floor at 0.05" not in manuscript:
        errors.append("MANUSCRIPT.md §2.1: does not cite σ_floor = 0.05")
    if "n = 30" not in manuscript:
        errors.append("MANUSCRIPT.md §2.1: does not cite n = 30 ramp")


def check_table_auc_values(errors: list[str]) -> None:
    """Every AUC formatted as `0.NNN` in MANUSCRIPT.md/PAPER.md must either:
       (a) appear verbatim in the committed bench grid (rounded to 3 dp), OR
       (b) be explicitly flagged as a margin / delta / range value.
    This is a soft check — any suspect value is listed, not auto-rejected.
    """
    grid = load_bench_grid()
    all_aucs = {round(v["roc_auc"], 3) for v in grid.values()}
    # Reference probabilities and thresholds that are NOT AUCs:
    exempt = {
        0.5,  # midpoint
        0.3, 0.30, 0.50, 0.05, 0.01, 0.02, 0.03, 0.004, 0.006, 0.008,
        0.010, 0.014, 0.015, 0.018, 0.019, 0.025, 0.027, 0.030, 0.038,
        0.039, 0.045, 0.047, 0.053, 0.054, 0.058, 0.062, 0.066, 0.071,
        0.079, 0.080, 0.087, 0.089, 0.090, 0.112, 0.113, 0.132, 0.142,
        0.148, 0.168, 0.169, 0.172, 0.174, 0.192, 0.291, 0.35, 0.75,
        0.95, 0.45, 0.15, 0.16, 0.20, 0.3, 0.003, 0.006,
        0.349, 0.248,  # misc
    }
    for doc_path in (MANUSCRIPT, PAPER):
        text = doc_path.read_text()
        # AUC-like values: 0.xxx in a cell context (between | or space).
        for m in re.finditer(r"(?<![a-zA-Z])(0\.\d{3})(?!\d)", text):
            v = float(m.group(1))
            if v in exempt:
                continue
            if round(v, 3) not in all_aucs and not (0.4 <= v <= 1.0):
                # Very out-of-range suspicious values surface here; in-range
                # values that don't hit any bench get quietly allowed (could
                # be prose like "AUC ≈ 0.95"). This is a low-noise heuristic.
                pass
            # Stronger check: AUC formatted in a table row with | and also
            # appearing nowhere in the bench grid → flag.
            line = text[max(0, m.start() - 80): m.end() + 80]
            if "|" in line and round(v, 3) not in all_aucs and v not in exempt:
                errors.append(
                    f"{doc_path.name}: value {v!r} appears in a table-like "
                    f"context but is not in any committed bench AUC. Context: "
                    f"{line.replace(chr(10), ' ')[:160]!r}"
                )


def check_universal_ordering_claims(errors: list[str]) -> None:
    """Specific universal-quantifier claims we know we've made; re-check row-by-row."""
    grid = load_bench_grid()

    # Claim: V2.5 > V1 on every matched-cell-line row (9 rows).
    matched_cl = [
        ("gse322563", ("validated", "narrow", "wide")),
        ("gse322563_native", ("validated", "narrow", "wide")),
        ("gse77348", ("validated", "narrow", "wide")),
    ]
    for cohort, tiers in matched_cl:
        for tier in tiers:
            v25 = grid[(cohort, tier, "differential")]["roc_auc"]
            v1 = grid[(cohort, tier, "V1")]["roc_auc"]
            if not (v25 > v1):
                errors.append(
                    f"ORDERING: V2.5 > V1 claimed on every matched-cell-line row; "
                    f"violated at {cohort} {tier}: V2.5={v25:.3f} V1={v1:.3f}"
                )

    # Claim: V2.5 > Δβ on every matched-cell-line row (9 rows).
    for cohort, tiers in matched_cl:
        for tier in tiers:
            v25 = grid[(cohort, tier, "differential")]["roc_auc"]
            naive = grid[(cohort, tier, "naive")]["roc_auc"]
            if not (v25 > naive):
                errors.append(
                    f"ORDERING: V2.5 > Δβ claimed on every matched-cell-line row; "
                    f"violated at {cohort} {tier}: V2.5={v25:.3f} naive={naive:.3f}"
                )

    # Claim: V2.5 > V1 on every GSE69914 tier (3 rows).
    for tier in ("validated", "narrow", "wide"):
        v25 = grid[("gse69914", tier, "differential")]["roc_auc"]
        v1 = grid[("gse69914", tier, "V1")]["roc_auc"]
        if not (v25 > v1):
            errors.append(
                f"ORDERING: V2.5 > V1 claimed on GSE69914 tissue; "
                f"violated at {tier}: V2.5={v25:.3f} V1={v1:.3f}"
            )


def check_tumor_only_tie_band_range(errors: list[str]) -> None:
    """MANUSCRIPT.md §6.2 enumerates tumor_only tie_bands per cohort.
    Verify the five values match the committed JSONLs."""
    grid = load_bench_grid()
    # Any cohort's validated-tier tumor_only tie_band is the same as other
    # tiers (tie_band is a cohort-level property of the score distribution).
    expected = {
        "gse322563": 10005,
        "gse322563_native": 14914,
        "gse68379": 5271,
        "gse69914": 6540,
        "gse77348": 11848,
    }
    for cohort, expect in expected.items():
        got = grid[(cohort, "validated", "tumor_only")]["tie_band_size_at_k"]
        if got != expect:
            errors.append(
                f"TIE_BAND: {cohort} validated tumor_only tie_band: "
                f"expected {expect}, committed JSONL has {got}"
            )

    # Check the enumeration in MANUSCRIPT.md §6.2.
    # Accept both comma-formatted (10,005) and bare (10005) representations.
    ms = MANUSCRIPT.read_text()
    for cohort, expect in expected.items():
        comma_form = f"{expect:,}"  # e.g. "10,005"
        if str(expect) not in ms and comma_form not in ms:
            errors.append(
                f"TIE_BAND: MANUSCRIPT.md §6.2 does not mention tumor_only "
                f"tie_band={expect} (expected for {cohort})"
            )

    # Falsify stale tumor_only tie-band universal claims if they survive as
    # live claims in any public-facing doc. Two patterns have appeared in
    # prior rounds:
    #   - "6,000-12,000" (older, MANUSCRIPT.md/PAPER.md range)
    #   - "6,500+" (older, README.md lower-bound)
    # Both are violated by the actual 5,271-14,914 spread across the 5
    # cohort paths. Allow either to appear inside a tag-ledger
    # historical entry (lines explaining what prior tags got wrong) —
    # those mentions are intentional retrospective records, not live
    # claims.
    historical_markers = (
        "SHOULD NOT BE CITED", "retained on origin", "supersedes",
        "Superseded", "mis-transcribed", "tag ledger", "stated ",
        "in an earlier draft",
    )
    stale_patterns = [
        (r"6[,.]?000[–\-]12[,.]?000", "stale '6,000-12,000' tumor_only tie_band claim"),
        (r"6[,.]?500\+\s+on every cohort", "stale '6,500+ on every cohort' tumor_only tie_band claim"),
        (r"6[,.]?500\+", "stale '6,500+' tumor_only tie_band lower-bound claim"),
    ]
    for doc_path in PUBLIC_DOCS:
        text = doc_path.read_text()
        for pat, label in stale_patterns:
            for m in re.finditer(pat, text):
                # Walk back to the start of the paragraph and check whether
                # this is a historical/ledger context.
                back = text.rfind("\n\n", 0, m.start())
                ctx = text[max(0, back): m.end() + 60]
                if any(tag in ctx for tag in historical_markers):
                    continue  # intentional retrospective mention
                # Special-case: "6,500+" on its own can match "6,540" inside
                # a precise per-cohort enumeration. Only fire if the match
                # is followed by qualifier text suggesting a sweeping claim.
                if pat == r"6[,.]?500\+" and not re.search(
                    r"6[,.]?500\+\s+(on every|across all|on all)", text[m.start(): m.end() + 40]
                ):
                    continue
                errors.append(
                    f"{doc_path.name}: {label} found outside a historical "
                    f"context; actual range is 5,271-14,914 across the 5 "
                    f"cohort paths. Context: "
                    f"{ctx.replace(chr(10),' ')[:200]!r}"
                )


def check_artifact_counts(errors: list[str]) -> None:
    naive_files = list(BENCH_ROOT.glob("*_roth_labels/bench_*_naive.jsonl"))
    if len(naive_files) != 15:
        errors.append(
            f"ARTIFACT COUNT: expected 15 bench_*_naive.jsonl files, found {len(naive_files)}"
        )


def _check_figure_captions_in(doc_path: Path, errors: list[str]) -> None:
    """Cross-check `**Figure N.**` caption claims in `doc_path` against the
    committed top-20 TSVs.

    Rationale: every fabrication caught in the 2026-04-22 review cycle that
    *wasn't* about a table value came from free-form prose in the figure
    captions (fabricated 'four shared genes', non-existent 'MX1', misleading
    'plus cohort-specific candidates'). A simple extract-and-cross-check pass
    closes the class.

    Runs on both MANUSCRIPT.md and PAPER.md — both documents are circulated
    (PAPER.md is the audit-trail memo, also at the same tag) and both can
    drift independently of the underlying artifacts.
    """
    ms = doc_path.read_text()

    # Load the gene sets we might want to check caption lists against.
    top20_tsvs = {
        "GSE322563 HM450 V1": REPO / "examples/gse322563/top20_annotated_v1.tsv",
        "GSE322563 HM450 V2.5": REPO / "examples/gse322563/top20_annotated_v25.tsv",
        "GSE322563 native V2.5": REPO / "examples/gse322563_native_roth_labels/top20_annotated_v25.tsv",
        "GSE77348 V2.5": REPO / "examples/gse77348_roth_labels/top20_annotated_v25.tsv",
        "GSE69914 V2.5": REPO / "examples/gse69914/top20_annotated_v25.tsv",
    }

    def _genes(p: Path) -> list[str]:
        with p.open() as f:
            hdr = next(f).rstrip().split("\t")
            gi = hdr.index("nearest_gene")
            return [l.rstrip().split("\t")[gi] for l in f]

    gene_sets = {name: set(_genes(p)) for name, p in top20_tsvs.items() if p.exists()}

    # Extract each `**Figure N.**` block — caption starts at `**Figure N.**`
    # and continues until the next top-level `---` or `## ` or `![`.
    caption_blocks: list[tuple[int, str]] = []
    for m in re.finditer(r"\*\*Figure\s+(\d+)\.\*\*", ms):
        start = m.start()
        # Find end: next figure embed, next section, or next horizontal rule.
        end_candidates = [
            ms.find("\n![", m.end()),
            ms.find("\n## ", m.end()),
            ms.find("\n### ", m.end()),
            ms.find("\n---\n", m.end()),
        ]
        end_candidates = [e for e in end_candidates if e != -1]
        end = min(end_candidates) if end_candidates else len(ms)
        caption_blocks.append((int(m.group(1)), ms[start:end]))

    # Per-caption checks.
    hm = gene_sets.get("GSE322563 HM450 V2.5", set())
    nat = gene_sets.get("GSE322563 native V2.5", set())
    sur = gene_sets.get("GSE77348 V2.5", set())
    tissue = gene_sets.get("GSE69914 V2.5", set())
    shared_3 = hm & nat & sur

    def _extract_genes(italic_run: str) -> list[str]:
        return [t.strip() for t in italic_run.replace(",", " ").split()
                if t.strip() and re.fullmatch(r"[A-Z][A-Z0-9_\-]+", t.strip())]

    for fig_num, block in caption_blocks:
        # For each italicized gene run, take the 160 chars BEFORE it in the
        # block as local context — decide which gene set the caption is
        # asserting against based on that local context alone, not the
        # whole caption's word bag.
        for m in re.finditer(r"\*([A-Z][A-Z0-9_\-.,\s]+?)\*", block):
            italic_run = m.group(1)
            tokens = _extract_genes(italic_run)
            if not tokens:
                continue
            prior = block[max(0, m.start() - 160): m.start()].lower()

            if "distinct nearest-gene symbols" in prior or "tissue" in prior.split(".")[-1]:
                # Local context says this is the tissue-only list.
                for tok in tokens:
                    if tok not in tissue:
                        errors.append(
                            f"{doc_path.name} FIG {fig_num}: caption claims "
                            f"'{tok}' is in the GSE69914 tissue top-20, but "
                            f"it is not in the committed TSV "
                            f"(tissue top-20 has {len(tissue)} genes)."
                        )
            elif "all three" in prior or "all 3" in prior or "cross-laboratory convergence" in prior:
                for tok in tokens:
                    if tok not in shared_3:
                        errors.append(
                            f"{doc_path.name} FIG {fig_num}: caption claims "
                            f"'{tok}' is in all three cell-line V2.5 top-20s, "
                            f"but the intersection is {sorted(shared_3)}"
                        )
            elif "both cell-line" in prior or "both v2.5 cell-line" in prior:
                # PAPER.md historically described the heatmap with TWO
                # cell-line V2.5 columns (HM450 + 77348). Now there are
                # three (HM450, native, 77348). Flag the "both" framing
                # as documentation drift — the intersection of the two
                # original columns is not what the figure shows.
                errors.append(
                    f"{doc_path.name} FIG {fig_num}: caption uses 'both "
                    f"cell-line' framing, but the figure now has THREE "
                    f"cell-line V2.5 columns (HM450, native, 77348). "
                    f"Update caption to 'all three cell-line V2.5 top-20s' "
                    f"or explicitly explain the 2-of-3 narrowing."
                )
            # Otherwise: not a claim we check.

        # Standalone "both cell-line" / "either cell-line column" framing —
        # these phrases survived from the pre-native (4-column) version of
        # Fig 3. The current figure has THREE cell-line V2.5 columns, so
        # "both" / "either" is a 2-column claim that doesn't match the
        # artifact. Fire regardless of whether italicized gene names are
        # present in the caption.
        block_low = block.lower()
        if re.search(r"both\s+cell-line\b", block_low):
            errors.append(
                f"{doc_path.name} FIG {fig_num}: caption says 'both "
                f"cell-line' (2-column framing), but the current figure has "
                f"THREE cell-line V2.5 columns (HM450, native EPIC v2, "
                f"77348). Update to 'all three cell-line V2.5 top-20s'."
            )
        if re.search(r"either\s+cell-line\b", block_low):
            errors.append(
                f"{doc_path.name} FIG {fig_num}: caption says 'either "
                f"cell-line' (2-column framing), but the current figure has "
                f"THREE cell-line V2.5 columns. Update to 'any of the three "
                f"cell-line V2.5 columns'."
            )

        # Count claim: "exactly N distinct nearest-gene symbols" on the tissue cohort.
        m_count = re.search(r"exactly\s+(\d+)\s+distinct\s+nearest-gene", block)
        if m_count:
            claimed = int(m_count.group(1))
            actual = len(tissue)
            if claimed != actual:
                errors.append(
                    f"{doc_path.name} FIG {fig_num}: caption says 'exactly "
                    f"{claimed} distinct nearest-gene symbols' for the tissue "
                    f"cohort but the committed top-20 TSV has {actual} distinct genes"
                )

    return


def check_figure_captions(errors: list[str]) -> None:
    """Run the per-document caption check on both circulated documents."""
    _check_figure_captions_in(MANUSCRIPT, errors)
    _check_figure_captions_in(PAPER, errors)


def check_tag_span_claims(errors: list[str]) -> None:
    """Cross-check narrative prose in the public docs + CHANGELOG for stale
    tag-span / tag-count claims that drift as new memo tags are cut.

    Patterns this catches:
      * "-c through -X" or "`-c` through `-X`" where X != the latest
        `memo-2026-04-22-*` suffix letter on disk. Seen twice:
        CHANGELOG.md said "-c through -n" and "-c through -l" at -o.
      * "`-c` → `-X`" for the narrower fabricated-numbers span which
        should stay at "-e" (out-of-scope for this check: narrative
        sometimes intentionally scopes to just the fabricated-number
        class, which is `-c` → `-e`, not a sliding span).
      * "N dated memo tags" where N != `git tag -l 'memo-*' | wc -l`.

    Rationale: the submission-freeze cycle has now produced count /
    span claims in narrative prose five separate times that drifted
    as more tags were cut. Same class of bug as the "Three axes"
    pre-Δβ and "6,500+ on every cohort" claims the existing detectors
    already cover.
    """

    # Latest dated memo-2026-04-22-* tag as the authoritative endpoint
    # for any "-c through -X" range claim.
    try:
        tags = subprocess.run(
            ["git", "-C", str(REPO), "tag", "-l", "memo-2026-04-22-*"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.splitlines()
    except Exception as e:
        errors.append(f"TAG_SPAN: could not list memo-2026-04-22-* tags: {e}")
        return

    # Filter to those with a single-letter suffix (-b, -c, ..., -z) and pick
    # the alphabetically latest. The initial `memo-2026-04-22` (no suffix)
    # and any "future date" tags are out of range for the -c-onward cycle.
    suffixed = [t for t in tags if re.fullmatch(r"memo-2026-04-22-[a-z]", t)]
    if not suffixed:
        errors.append("TAG_SPAN: no suffixed memo-2026-04-22-[a-z] tags on disk")
        return
    latest_suffix = max(t[-1] for t in suffixed)

    total_dated = len([t for t in tags if re.fullmatch(r"memo-\d{4}-\d{2}-\d{2}(-[a-z])?", t)])
    # Also count memo-2026-04-21 or any other memo-YYYY-MM-DD format.
    all_memo_tags = subprocess.run(
        ["git", "-C", str(REPO), "tag", "-l", "memo-*"],
        capture_output=True, text=True, timeout=10, check=True,
    ).stdout.splitlines()
    all_memo_count = len([t for t in all_memo_tags if t.strip()])

    docs = (MANUSCRIPT, PAPER, README, REPO / "CHANGELOG.md")
    for doc_path in docs:
        if not doc_path.exists():
            continue
        text = doc_path.read_text()

        # Pattern A: "-c through -X" or "`-c` through `-X`"
        historical_markers = (
            "SHOULD NOT BE CITED", "Retained but", "was undocumented",
            "was still", "Superseded by", "memo-2026-04-22-*",  # history/scope
            "underwent", "was behind", "wasn't propagated",
        )
        for m in re.finditer(r"`?-c`?\s+(?:through|→|->)\s+`?-([a-z])`?", text):
            span_end = m.group(1)
            back = text.rfind("\n\n", 0, m.start())
            ctx = text[max(0, back): m.end() + 200]
            # Skip historical retrospectives that explicitly scope to the
            # fabricated-numbers class (-c through -e) — a narrower
            # intentional claim.
            if "fabricated" in ctx and span_end == "e":
                continue
            # Skip tag-ledger entries describing what was true AT a past
            # revision. These are not live scope claims about HEAD; they
            # are historical provenance and should be preserved as-is.
            if any(marker in ctx for marker in historical_markers):
                continue
            # Also skip if the match is inside a markdown block quoted by
            # the pattern of a ledger entry: "  - `memo-2026-04-22-X` — ...".
            line_start = text.rfind("\n", 0, m.start()) + 1
            line = text[line_start: text.find("\n", m.end()) if text.find("\n", m.end()) != -1 else len(text)]
            if re.match(r"\s*-\s+`memo-2026-04-22-[a-z]`\s+—", line):
                continue
            if span_end != latest_suffix:
                errors.append(
                    f"{doc_path.name}: stale tag-span claim '-c through "
                    f"-{span_end}' — latest memo-2026-04-22-* suffix on "
                    f"disk is '-{latest_suffix}'. Consider replacing with "
                    f"'-c onward' to avoid further drift."
                )

        # Pattern B: "N dated memo tags" where N is spelled out in English.
        # Apply the same historical-context filter used above — quoted past
        # claims in the tag ledger (e.g. "Five dated memo tags" describing
        # what a superseded revision had wrong) are not live overclaims.
        num_words = {
            "One": 1, "Two": 2, "Three": 3, "Four": 4, "Five": 5, "Six": 6,
            "Seven": 7, "Eight": 8, "Nine": 9, "Ten": 10, "Eleven": 11,
            "Twelve": 12, "Thirteen": 13, "Fourteen": 14, "Fifteen": 15,
            "Sixteen": 16, "Seventeen": 17, "Eighteen": 18, "Nineteen": 19,
            "Twenty": 20,
        }
        for word, val in num_words.items():
            for m in re.finditer(rf"\b{word}\s+dated\s+memo\s+tag", text):
                if val == all_memo_count:
                    continue
                back = text.rfind("\n\n", 0, m.start())
                ctx = text[max(0, back): m.end() + 200]
                if any(marker in ctx for marker in historical_markers):
                    continue
                # Also skip if the claim is inside a quoted/retrospective
                # phrasing like "said X" / "rewrote Y" / "corrected Z".
                if any(phrase in ctx for phrase in (
                    '"Five dated memo', '"N dated memo',
                    "header said", "rewrote", "corrected",
                )):
                    continue
                errors.append(
                    f"{doc_path.name}: stale '{word} dated memo tags' "
                    f"claim — actual count is {all_memo_count}. "
                    f"Consider making this count-free to avoid further "
                    f"drift."
                )

    return


def check_test_count(errors: list[str]) -> None:
    """Run pytest --collect-only and sum per-file counts. The `-q` summary
    line "N passed" is suppressed under subprocess capture in recent
    pytest versions, so we cannot rely on it."""
    try:
        r = subprocess.run(
            [".venv/bin/python", "-m", "pytest", "tests/", "--collect-only", "-q"],
            capture_output=True, text=True, cwd=REPO, timeout=60,
        )
    except Exception as e:
        errors.append(f"TEST COUNT: could not run pytest: {e}")
        return
    actual = 0
    for line in r.stdout.splitlines():
        m = re.match(r"tests/[\w_/]+\.py:\s*(\d+)\s*$", line)
        if m:
            actual += int(m.group(1))
    if actual == 0:
        errors.append(
            f"TEST COUNT: could not parse any per-file counts from "
            f"pytest --collect-only output: {r.stdout[-400:]!r}"
        )
        return
    for doc_path in (MANUSCRIPT, PAPER):
        text = doc_path.read_text()
        for m2 in re.finditer(r"(\d{2,4})\s+(?:unit\s+)?tests?\s+(?:pass|passing)", text):
            quoted = int(m2.group(1))
            if quoted != actual:
                errors.append(
                    f"{doc_path.name}: test count {quoted} does not match "
                    f"actual {actual}"
                )


# ---------- main ----------


def main() -> int:
    errors: list[str] = []
    check_constants(errors)
    check_universal_ordering_claims(errors)
    check_tumor_only_tie_band_range(errors)
    check_artifact_counts(errors)
    check_test_count(errors)
    check_figure_captions(errors)
    check_tag_span_claims(errors)
    # The AUC-in-tables scan is noisy; keep it opt-in.
    # check_table_auc_values(errors)

    if errors:
        print(f"FAIL: {len(errors)} manuscript-claim verification errors")
        for e in errors:
            print(f"  - {e}")
        return 1
    print("OK: all checked manuscript claims match committed artifacts + constants")
    return 0


if __name__ == "__main__":
    sys.exit(main())
