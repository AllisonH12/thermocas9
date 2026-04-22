#!/usr/bin/env bash
#
# Render PAPER.md → PAPER.pdf as a preprint.
#
# Toolchain (install once):
#     brew install pandoc typst
#
# Pandoc converts the markdown to typst syntax; typst compiles to PDF.
# The combined install is ~315 MB and avoids a full LaTeX distribution.
#
# Run from the repo root:
#     bash scripts/render_paper_pdf.sh
#
# Output: PAPER.pdf at the repo root.
#
# What this script does on top of a vanilla `pandoc PAPER.md -o PAPER.pdf`:
#
#   * Strips the leading title block (H1 + author/date/code/status
#     paragraph) from a temp copy of PAPER.md so that pandoc does NOT
#     treat the title as section 1, which would push every real section
#     down one level and double-number them ("1.2.1 1.1 ThermoCas9...").
#   * Passes title / author / date / subtitle to pandoc as metadata so
#     the PDF gets a real title block.
#   * Adds a numbered table of contents.
#   * Lets the manual "1 · Background", "2 · The V2 misspecification",
#     etc. headings carry their own numbering (no pandoc auto-numbering),
#     so the TOC matches what the markdown body says.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${REPO_ROOT}/PAPER.md"
OUT="${REPO_ROOT}/PAPER.pdf"

if [ ! -f "$SRC" ]; then
    echo "PAPER.md not found at $SRC" >&2
    exit 1
fi
for tool in pandoc typst; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not on PATH; install with: brew install pandoc typst" >&2
        exit 1
    fi
done

TMP="$(mktemp -t paper.XXXXXX.md)"
trap 'rm -f "$TMP"' EXIT

# Strip the leading title block:
#   * delete the first H1 line ("# ...")
#   * delete the immediately-following author/date/code/status paragraph
#     (lines starting with **Author.** through the next blank line)
# Everything else is preserved verbatim.
awk '
  BEGIN { in_meta = 0; stripped_h1 = 0 }
  # First H1 → drop, enter meta-strip state
  !stripped_h1 && /^# / { stripped_h1 = 1; in_meta = 1; next }
  # While in meta-strip state, drop lines until the first blank line
  in_meta {
      if ($0 == "") { in_meta = 0 }
      next
  }
  { print }
' "$SRC" > "$TMP"

pandoc "$TMP" -o "$OUT" --pdf-engine=typst \
    --metadata title="Differential-protection probabilistic scoring for methylome-guided ThermoCas9 target-site ranking" \
    --metadata subtitle="Technical memo from the thermocas methylome-analysis framework" \
    --metadata author="Allison Huang (Thermocas9 Inc)" \
    --metadata date="$(date +%Y-%m-%d)" \
    --toc --toc-depth=2 \
    --variable=fontsize=10pt

ls -lh "$OUT"
