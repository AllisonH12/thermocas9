#!/usr/bin/env bash
#
# Render a markdown source document to PDF as a preprint.
#
# Usage:
#     bash scripts/render_paper_pdf.sh                    # default: PAPER.md  → PAPER.pdf
#     bash scripts/render_paper_pdf.sh PAPER.md           # explicit
#     bash scripts/render_paper_pdf.sh MANUSCRIPT.md      # → MANUSCRIPT.pdf
#     bash scripts/render_paper_pdf.sh path/to/doc.md     # → path/to/doc.pdf
#
# The submission-shaped Bioinformatics manuscript is `MANUSCRIPT.md`;
# the longer audit-trail memo is `PAPER.md`. They share the same
# rendering pipeline but produce distinct PDFs and metadata blocks.
#
# Toolchain (install once):
#     brew install pandoc typst
#
# Pandoc converts the markdown to typst syntax; typst compiles to PDF.
# The combined install is ~315 MB and avoids a full LaTeX distribution.
#
# What this script does on top of `pandoc <src>.md -o <src>.pdf`:
#
#   * Strips the leading title block (H1 + author/date/code/status
#     paragraph) from a temp copy so that pandoc does NOT treat the
#     title as section 1, which would push every real section down
#     one level and double-number them.
#   * Passes title / author / date / subtitle to pandoc as metadata
#     so the PDF gets a real title block.
#   * Adds a numbered table of contents.
#   * Lets the manual section headings carry their own numbering.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC_REL="${1:-PAPER.md}"

# Resolve to absolute. Allow the caller to pass either a path relative
# to the repo root (e.g. "MANUSCRIPT.md") or an absolute path.
if [ -f "$SRC_REL" ]; then
    SRC="$(cd "$(dirname "$SRC_REL")" && pwd)/$(basename "$SRC_REL")"
elif [ -f "${REPO_ROOT}/${SRC_REL}" ]; then
    SRC="${REPO_ROOT}/${SRC_REL}"
else
    echo "source markdown not found: $SRC_REL (looked in CWD and ${REPO_ROOT})" >&2
    exit 1
fi

# Output: same dir + same stem as source, with .pdf extension.
SRC_DIR="$(dirname "$SRC")"
SRC_STEM="$(basename "$SRC" .md)"
OUT="${SRC_DIR}/${SRC_STEM}.pdf"

for tool in pandoc typst; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not on PATH; install with: brew install pandoc typst" >&2
        exit 1
    fi
done

# Per-document metadata. Both share the author / date; title and subtitle
# differ. Add new doc stems to this case as needed.
case "$SRC_STEM" in
    PAPER)
        TITLE="Differential-protection probabilistic scoring for methylome-guided ThermoCas9 target-site ranking"
        SUBTITLE="Technical memo from the thermocas methylome-analysis framework (audit-trail revision)"
        ;;
    MANUSCRIPT)
        TITLE="thermocas: probabilistic ranking of methylation-protected ThermoCas9 target sites with tie-aware benchmarking"
        SUBTITLE="Bioinformatics-submission-shaped manuscript"
        ;;
    *)
        TITLE="$SRC_STEM"
        SUBTITLE=""
        ;;
esac
AUTHOR="Allison Huang (Columbia University)"

# mktemp -t on macOS doesn't preserve the trailing .md, so pandoc can't
# deduce the format from extension. Create a directory and put a .md
# file inside it — portable and avoids the warning.
TMPDIR="$(mktemp -d -t pdf_render.XXXXXX)"
TMP="${TMPDIR}/src.md"
trap 'rm -rf "$TMPDIR"' EXIT

# Strip the leading title block:
#   * delete the first H1 line ("# ...")
#   * delete the immediately-following author/date/code/status paragraph
#     (lines until the first blank line)
# Everything else is preserved verbatim.
awk '
  BEGIN { in_meta = 0; stripped_h1 = 0 }
  !stripped_h1 && /^# / { stripped_h1 = 1; in_meta = 1; next }
  in_meta {
      if ($0 == "") { in_meta = 0 }
      next
  }
  { print }
' "$SRC" > "$TMP"

# Resource path: figures referenced via relative paths
# (e.g. ![Figure 1](docs/figures/fig1...)) need to resolve from the
# repo root, not from /tmp where TMP lives. --resource-path tells
# pandoc where to look for those.
pandoc "$TMP" -o "$OUT" --pdf-engine=typst \
    --resource-path="${REPO_ROOT}" \
    --metadata title="$TITLE" \
    --metadata subtitle="$SUBTITLE" \
    --metadata author="$AUTHOR" \
    --metadata date="$(date +%Y-%m-%d)" \
    --toc --toc-depth=2 \
    --variable=fontsize=10pt

ls -lh "$OUT"
