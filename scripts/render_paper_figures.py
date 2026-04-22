"""Render the three figures referenced in PAPER.md.

This is the only file in the repo with a hard dependency on matplotlib.
The package itself is stdlib-only; install matplotlib separately for
figure rendering:

    uv pip install matplotlib

Run from the repo root:

    python scripts/render_paper_figures.py

Outputs (PNG and SVG side-by-side):

    docs/figures/fig1_mode_schematic.{png,svg}
    docs/figures/fig2_auc_bars.{png,svg}
    docs/figures/fig3_topgene_heatmap.{png,svg}

Numbers are hard-coded from the committed BenchmarkResult JSONLs and
the committed top-20 annotated TSVs (examples/*) so the figures match
the values in PAPER.md exactly. If a benchmark is rerun and the
artifacts change, regenerate by re-running this script.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


# ---------- Figure 1 — V2.5 mode-formula schematic ----------


def fig1_mode_schematic() -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")

    # Title
    ax.text(5, 4.6, "V2.5 differential-protection probabilistic composite",
            ha="center", fontsize=13, fontweight="bold")

    # Three factor boxes (top row)
    box_color = "#e8f0fe"
    edge_color = "#1a73e8"
    factors = [
        (1.5, "p_targ",  "P(β_tumor < 0.30)\n(method-of-moments Beta)"),
        (5.0, "p_diff",  "P(β_normal − β_tumor > δ)\n(normal approximation)"),
        (8.5, "p_trust", "EvidenceClass × min(1, n/30)\n(saturating)"),
    ]
    for x, label, sub in factors:
        ax.add_patch(mpatches.FancyBboxPatch(
            (x - 1.1, 2.6), 2.2, 1.0,
            boxstyle="round,pad=0.08", linewidth=1.5,
            facecolor=box_color, edgecolor=edge_color))
        ax.text(x, 3.3, label, ha="center", va="center",
                fontsize=12, fontweight="bold", family="monospace")
        ax.text(x, 2.85, sub, ha="center", va="center", fontsize=8.5)

    # Multiplication signs
    for x in (3.05, 6.55):
        ax.text(x, 3.1, "×", ha="center", va="center",
                fontsize=20, color=edge_color, fontweight="bold")

    # Composite arrow + result
    ax.annotate("", xy=(5, 1.4), xytext=(5, 2.55),
                arrowprops=dict(arrowstyle="->", color=edge_color, lw=2))
    ax.add_patch(mpatches.FancyBboxPatch(
        (3.0, 0.55), 4.0, 0.85,
        boxstyle="round,pad=0.08", linewidth=1.5,
        facecolor="#fff5e6", edgecolor="#f9ab00"))
    ax.text(5, 1.10, "p_therapeutic_selectivity", ha="center", va="center",
            fontsize=11, fontweight="bold", family="monospace")
    ax.text(5, 0.78, "  =  p_targ × p_diff × p_trust", ha="center", va="center",
            fontsize=9.5, family="monospace")

    # Legend / margin note
    ax.text(0.05, 0.05,
            "δ = differential margin (default 0.2). σ_Δ = √(σ_t² + σ_n²) "
            "with σ_k = IQR_k / 1.349 (floor 0.05).",
            fontsize=8, color="#5f6368")

    plt.tight_layout()
    fig.savefig(OUT / "fig1_mode_schematic.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig1_mode_schematic.svg",       bbox_inches="tight")
    plt.close(fig)


# ---------- Figure 2 — AUC bar charts ----------


def _auc_from_bench(path: Path) -> float:
    return json.loads(path.read_text().strip())["roc_auc"]


def fig2_auc_bars() -> None:
    bench = REPO / "examples"
    cohorts = [
        ("GSE322563\n(Roth, n=2/2)",   "gse322563_roth_labels"),
        ("GSE77348\n(surrogate, n=3/3)", "gse77348_roth_labels"),
        ("GSE69914\n(tissue, n=305/50)",  "gse69914_roth_labels"),
    ]
    label_sets = ["validated", "narrow", "wide"]
    axes_modes = [
        ("V1 final_score",        "V1",          "#5f6368"),
        ("V2 tumor_only",         "tumor_only",  "#fbbc04"),
        ("V2.5 differential",     "differential", "#1a73e8"),
    ]

    auc = {}
    for cohort_label, dirname in cohorts:
        for ls in label_sets:
            for axis_label, mode_tag, _ in axes_modes:
                p = bench / dirname / f"bench_{ls}_{mode_tag}.jsonl"
                auc[(cohort_label, ls, axis_label)] = _auc_from_bench(p)

    # Two-panel layout: left = primary endpoint (validated only), right = sensitivity
    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(12.5, 5.0), gridspec_kw={"width_ratios": [1, 2]}
    )

    # Left panel: validated AUC, 3 cohorts × 3 axes (9 bars in 3 groups)
    width = 0.27
    x = list(range(len(cohorts)))
    for i, (axis_label, _, color) in enumerate(axes_modes):
        ys = [auc[(c[0], "validated", axis_label)] for c in cohorts]
        offsets = [xi + (i - 1) * width for xi in x]
        axL.bar(offsets, ys, width=width, color=color, label=axis_label,
                edgecolor="#202124", linewidth=0.4)
        for off, y in zip(offsets, ys):
            axL.text(off, y + 0.012, f"{y:.2f}", ha="center", fontsize=7.5)
    axL.axhline(0.5, color="#dadce0", linestyle="--", linewidth=1)
    axL.set_xticks(x)
    axL.set_xticklabels([c[0] for c in cohorts], fontsize=8.5)
    axL.set_ylabel("AUC")
    axL.set_ylim(0, 1.05)
    axL.set_title("(a) Primary + tissue: validated label set (n_pos=3)\n"
                  "AUC at Roth Fig. 5d target probes",
                  fontsize=10)
    axL.legend(loc="lower left", fontsize=8, frameon=False)
    axL.spines["top"].set_visible(False)
    axL.spines["right"].set_visible(False)

    # Right panel: sensitivity over label granularity (3 cohorts × 3 sets × 3 axes)
    # Group by (cohort, label_set), 3 bars per group.
    groups = [(c[0], ls) for c in cohorts for ls in label_sets]
    x2 = list(range(len(groups)))
    for i, (axis_label, _, color) in enumerate(axes_modes):
        ys = [auc[(c, ls, axis_label)] for (c, ls) in groups]
        offsets = [xi + (i - 1) * width for xi in x2]
        axR.bar(offsets, ys, width=width, color=color,
                edgecolor="#202124", linewidth=0.4)
    axR.axhline(0.5, color="#dadce0", linestyle="--", linewidth=1)
    axR.set_xticks(x2)
    axR.set_xticklabels([f"{c.split(chr(10))[0]}\n{ls}" for (c, ls) in groups],
                        fontsize=7.0)
    axR.set_ylim(0, 1.05)
    axR.set_title("(b) Sensitivity: AUC across label granularities\n"
                  "(narrow ±50 bp, wide ±500 bp)",
                  fontsize=10)
    axR.spines["top"].set_visible(False)
    axR.spines["right"].set_visible(False)

    fig.suptitle("Figure 2 · Cross-cohort AUC by score axis (primary cohorts)",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / "fig2_auc_bars.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig2_auc_bars.svg",       bbox_inches="tight")
    plt.close(fig)


# ---------- Figure 3 — top-20 gene presence heatmap ----------


def _genes_from_tsv(path: Path) -> list[str]:
    """Return ordered gene symbols (rank 1 first) from a top20_annotated TSV."""
    out = []
    with path.open() as f:
        header = next(f).rstrip("\n").split("\t")
        gi = header.index("nearest_gene")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            out.append(parts[gi])
    return out


def fig3_topgene_heatmap() -> None:
    examples = REPO / "examples"
    # (header_label, path, tie_band, header_color)
    columns = [
        ("GSE322563\nV1\ntie_band=1",        examples / "gse322563" / "top20_annotated_v1.tsv",            1, "#1a73e8"),
        ("GSE322563\nV2.5\ntie_band=190",    examples / "gse322563" / "top20_annotated_v25.tsv",         190, "#d93025"),
        ("GSE77348\nV2.5\ntie_band=299",     examples / "gse77348_roth_labels" / "top20_annotated_v25.tsv", 299, "#d93025"),
        ("GSE69914\nV2.5\ntie_band=2",       examples / "gse69914" / "top20_annotated_v25.tsv",            2, "#1a73e8"),
    ]
    col_genes = [(_genes_from_tsv(p), tb, label, color) for (label, p, tb, color) in columns]

    # Union of unique genes; sort so that most-shared rows appear at the top.
    all_genes_ordered: list[str] = []
    seen: set[str] = set()
    for genes, _, _, _ in col_genes:
        for g in genes:
            if g and g not in seen:
                seen.add(g)
                all_genes_ordered.append(g)

    presence = {
        g: sum(1 for genes, _, _, _ in col_genes if g in genes)
        for g in all_genes_ordered
    }
    first_col = {
        g: next(i for i, (genes, _, _, _) in enumerate(col_genes) if g in genes)
        for g in all_genes_ordered
    }
    rows = sorted(all_genes_ordered, key=lambda g: (-presence[g], first_col[g], g))
    matrix = [[1 if g in col_genes[c][0] else 0 for c in range(len(col_genes))] for g in rows]
    shared_cell_line_rows = [
        i for i, g in enumerate(rows)
        if g in col_genes[1][0] and g in col_genes[2][0]
    ]

    # Generous height — 0.22 in per gene row + headroom for the column header.
    fig_h = max(8.5, len(rows) * 0.22 + 2.0)
    fig, ax = plt.subplots(figsize=(8.5, fig_h))

    cmap = ListedColormap(["#f1f3f4", "#1a73e8"])
    ax.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest")

    ax.set_xticks(range(len(col_genes)))
    ax.set_xticklabels([])  # use ax.text to color column headers individually
    ax.tick_params(axis="x", length=0)
    for c, (_, _, label, color) in enumerate(col_genes):
        ax.text(c, -1.6, label, ha="center", va="center",
                fontsize=8.5, color=color, fontweight="bold",
                family="DejaVu Sans")

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=7.8)
    for i in shared_cell_line_rows:
        ax.get_yticklabels()[i].set_fontweight("bold")
        ax.get_yticklabels()[i].set_color("#1a73e8")

    # Tighter cell separators (gridlines between rows/columns).
    ax.set_xticks([x - 0.5 for x in range(1, len(col_genes))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(rows))], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)
    ax.set_xlim(-0.5, len(col_genes) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -3.2)

    # Title sits above the column headers, well clear of the heatmap body.
    fig.suptitle(
        "Figure 3 · Top-20 gene presence per (axis × cohort) — tie-window-aware",
        fontsize=11, fontweight="bold", y=0.985,
    )

    # Caveat block below the figure body.
    fig.text(
        0.02, -0.02,
        "Blue column headers: tie_band ≤ 2 — top-20 is the genuine top-20 of the score distribution.\n"
        "Red column headers:  tie_band ≫ K — top-20 is a 20-record window inside the tied band\n"
        "                      selected by the documented `candidate_id` ascending tie-break.\n"
        "Gene rows highlighted in bold blue appear in BOTH cell-line V2.5 top-20 windows\n"
        "(GSE322563, GSE77348). On those cohorts shared membership reflects window convergence,\n"
        "not robust ranking convergence; AUC (Fig. 2) is the stable claim at low n.",
        fontsize=7.5, color="#5f6368", family="DejaVu Sans",
    )

    plt.tight_layout()
    fig.savefig(OUT / "fig3_topgene_heatmap.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3_topgene_heatmap.svg",       bbox_inches="tight")
    plt.close(fig)


# ---------- main ----------


def main() -> int:
    print("rendering Fig 1 (mode schematic)...")
    fig1_mode_schematic()
    print("rendering Fig 2 (AUC bars)...")
    fig2_auc_bars()
    print("rendering Fig 3 (top-20 heatmap)...")
    fig3_topgene_heatmap()
    print(f"→ {OUT}/")
    for p in sorted(OUT.iterdir()):
        print(f"  {p.relative_to(REPO)}  ({p.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
