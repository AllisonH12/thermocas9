#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib>=3.8"]
# ///
"""Atlas dot-plot — per-positive whole-genome rank percentiles for the website.

Produces the headline website figure for the /atlas page:
  - x-axis: cohort (4 paths, ordered cell-line → tissue)
  - y-axis: WG rank percentile (top of the universe = 100)
  - one row per Roth Fig. 5d positive (ESR1 / EGFLAM / GATA3), color-coded
  - two markers per (cohort, gene): V2.5-diff (open circle) and V2.5-sigmoid
    (filled), connected by a thin line so the rank-lift direction is visible

The data is sourced verbatim from the existing materialized tables:
  - examples/genome_wide_panel.md   (V2.5-diff and V2.5-sigmoid WG ranks
    for all 4 cohorts × 3 positives)
  - examples/tissue_per_positive_wg_ranks.tsv  (additionally Δβ-only and
    limma-style for GSE69914 tissue, used for the supplementary tissue
    panel)

WG denominators are hard-coded from the panel SHA256 records cited in
PAPER.md §5.2.2 / MANUSCRIPT.md §1 abstract (HM450 19.8M; EPIC v2 35.4M;
GSE69914 uses HM450; GSE77348 uses HM450). Percentile = 1 - rank/N.

Output:
  docs/figures/fig4_atlas_per_positive_wg_percentile.{png,svg}
  docs/website/atlas/per_positive_wg_percentile.json   (structured data,
    for an interactive Observable Plot / D3 component to consume later)

Run:
  uv run python scripts/build_atlas_dotplot.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# WG denominators from PAPER.md §5.2.2.
N_HM450 = 19_787_820
N_EPIC_V2 = 35_400_000  # rounded to the abstract's "35.4M"; rank-percentile
                        # error from the rounding is < 1 part in 10⁵.

# (cohort_label, axis, gene, wg_rank, denominator)
# Source: examples/genome_wide_panel.md "Per-positive WG ranks" table.
ROWS = [
    # GSE322563 HM450
    ("GSE322563 HM450",     "V2.5-diff",    "ESR1",     17_661, N_HM450),
    ("GSE322563 HM450",     "V2.5-sigmoid", "ESR1",     12_845, N_HM450),
    ("GSE322563 HM450",     "V2.5-diff",    "EGFLAM",  460_160, N_HM450),
    ("GSE322563 HM450",     "V2.5-sigmoid", "EGFLAM",  488_531, N_HM450),
    ("GSE322563 HM450",     "V2.5-diff",    "GATA3",   180_297, N_HM450),
    ("GSE322563 HM450",     "V2.5-sigmoid", "GATA3",   185_420, N_HM450),
    # GSE322563 native EPIC v2
    ("GSE322563 EPIC v2",   "V2.5-diff",    "ESR1",      5_746, N_EPIC_V2),
    ("GSE322563 EPIC v2",   "V2.5-sigmoid", "ESR1",      4_037, N_EPIC_V2),
    ("GSE322563 EPIC v2",   "V2.5-diff",    "EGFLAM",  158_502, N_EPIC_V2),
    ("GSE322563 EPIC v2",   "V2.5-sigmoid", "EGFLAM",  167_809, N_EPIC_V2),
    ("GSE322563 EPIC v2",   "V2.5-diff",    "GATA3",    50_096, N_EPIC_V2),
    ("GSE322563 EPIC v2",   "V2.5-sigmoid", "GATA3",    51_498, N_EPIC_V2),
    # GSE77348
    ("GSE77348",            "V2.5-diff",    "ESR1",     24_935, N_HM450),
    ("GSE77348",            "V2.5-sigmoid", "ESR1",     20_154, N_HM450),
    ("GSE77348",            "V2.5-diff",    "EGFLAM",  855_023, N_HM450),
    ("GSE77348",            "V2.5-sigmoid", "EGFLAM",  904_130, N_HM450),
    ("GSE77348",            "V2.5-diff",    "GATA3",   184_710, N_HM450),
    ("GSE77348",            "V2.5-sigmoid", "GATA3",   183_854, N_HM450),
    # GSE69914 tissue
    ("GSE69914 tissue",     "V2.5-diff",    "ESR1",    514_603, N_HM450),
    ("GSE69914 tissue",     "V2.5-sigmoid", "ESR1",  2_200_263, N_HM450),
    ("GSE69914 tissue",     "V2.5-diff",    "EGFLAM",9_183_566, N_HM450),
    ("GSE69914 tissue",     "V2.5-sigmoid", "EGFLAM",4_944_861, N_HM450),
    ("GSE69914 tissue",     "V2.5-diff",    "GATA3", 3_490_723, N_HM450),
    ("GSE69914 tissue",     "V2.5-sigmoid", "GATA3", 1_026_284, N_HM450),
]

COHORTS = [
    "GSE322563 HM450",
    "GSE322563 EPIC v2",
    "GSE77348",
    "GSE69914 tissue",
]

GENES = ["ESR1", "EGFLAM", "GATA3"]
GENE_COLORS = {"ESR1": "#1b9e77", "EGFLAM": "#d95f02", "GATA3": "#7570b3"}


def to_percentile(rank: int, denom: int) -> float:
    return 100.0 * (1.0 - rank / denom)


def by_key(cohort: str, axis: str, gene: str) -> float | None:
    for c, a, g, r, n in ROWS:
        if (c, a, g) == (cohort, axis, gene):
            return to_percentile(r, n)
    return None


def render_dotplot(out_png: Path, out_svg: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=180)

    # Within each cohort, jitter the three genes horizontally so dumbbells
    # don't overlap. ±0.18 spread.
    cohort_x = {c: i for i, c in enumerate(COHORTS)}
    gene_offset = {"ESR1": -0.18, "EGFLAM": 0.0, "GATA3": +0.18}

    for cohort in COHORTS:
        for gene in GENES:
            x = cohort_x[cohort] + gene_offset[gene]
            p_diff = by_key(cohort, "V2.5-diff", gene)
            p_sig = by_key(cohort, "V2.5-sigmoid", gene)
            color = GENE_COLORS[gene]

            ax.plot([x, x], [p_diff, p_sig], color=color, linewidth=1.0,
                    alpha=0.55, zorder=1)
            ax.scatter([x], [p_diff], s=58, facecolors="white",
                       edgecolors=color, linewidths=1.6, zorder=3)
            ax.scatter([x], [p_sig], s=58, facecolors=color,
                       edgecolors=color, linewidths=1.6, zorder=3)

    ax.set_xticks(range(len(COHORTS)))
    ax.set_xticklabels(COHORTS, fontsize=10)
    ax.set_ylabel("Whole-genome rank percentile  (top of universe = 100)",
                  fontsize=10)
    ax.set_ylim(0, 102)
    # Drop the 99 / 100 ticks — they collide visually at the top of the
    # plot where the cell-line dumbbells stack. The 95 dotted reference
    # line still tells the reader where "top 5%" sits.
    ax.set_yticks([0, 25, 50, 75, 90, 95])
    ax.axhline(99, color="#999999", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.axhline(95, color="#999999", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.axhline(75, color="#cccccc", linewidth=0.5, linestyle=":", alpha=0.6)
    ax.set_title(
        "Per-positive whole-genome rank percentiles — Roth Fig. 5d positives "
        "across the four evaluated cohort paths",
        fontsize=11, pad=14,
    )
    ax.grid(axis="y", color="#eeeeee", linewidth=0.5, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Two-part legend: gene color + axis marker style.
    gene_handles = [
        mpatches.Patch(color=GENE_COLORS[g], label=g) for g in GENES
    ]
    axis_handles = [
        plt.Line2D([0], [0], marker="o", linestyle="None",
                   markerfacecolor="white", markeredgecolor="#444",
                   markersize=8, markeredgewidth=1.6,
                   label="V2.5-diff (open)"),
        plt.Line2D([0], [0], marker="o", linestyle="None",
                   markerfacecolor="#444", markeredgecolor="#444",
                   markersize=8, markeredgewidth=1.6,
                   label="V2.5-sigmoid (filled)"),
    ]
    leg1 = ax.legend(handles=gene_handles, loc="lower left",
                     frameon=False, fontsize=9, title="gene")
    ax.add_artist(leg1)
    ax.legend(handles=axis_handles, loc="lower right",
              frameon=False, fontsize=9, title="scoring axis")

    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def emit_json(out_json: Path) -> None:
    payload = {
        "schema_version": 1,
        "description": (
            "Per-positive whole-genome rank percentiles for the three Roth "
            "Fig. 5d positives across the four evaluated cohort paths, "
            "scored under V2.5-diff and V2.5-sigmoid. Source: "
            "examples/genome_wide_panel.md (PAPER.md §5.2.2)."
        ),
        "wg_denominators": {
            "HM450": N_HM450,
            "EPIC_v2": N_EPIC_V2,
        },
        "rows": [
            {
                "cohort": c,
                "axis": a,
                "gene": g,
                "wg_rank": r,
                "wg_denominator": n,
                "wg_rank_percentile": round(to_percentile(r, n), 4),
            }
            for c, a, g, r, n in ROWS
        ],
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    fig_dir = repo_root / "docs" / "figures"
    web_dir = repo_root / "docs" / "website" / "atlas"
    fig_dir.mkdir(parents=True, exist_ok=True)
    web_dir.mkdir(parents=True, exist_ok=True)

    png = fig_dir / "fig4_atlas_per_positive_wg_percentile.png"
    svg = fig_dir / "fig4_atlas_per_positive_wg_percentile.svg"
    js = web_dir / "per_positive_wg_percentile.json"

    render_dotplot(png, svg)
    emit_json(js)

    print(f"wrote {png.relative_to(repo_root)}")
    print(f"wrote {svg.relative_to(repo_root)}")
    print(f"wrote {js.relative_to(repo_root)}")


if __name__ == "__main__":
    main()
