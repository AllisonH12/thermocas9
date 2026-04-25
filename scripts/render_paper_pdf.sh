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
        TITLE="Compositional probability-scale scoring and tie-band-aware benchmarking for methylome-guided ThermoCas9 target-site ranking"
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

# Date sourcing — reproducible from the tagged source.
#   1. Look for "**Date.** YYYY-MM-DD" in the source MD title block.
#   2. Fall back to the most-recent local memo-YYYY-MM-DD[*] tag's date.
#   3. Last resort: wall-clock date with a stderr warning. Re-rendering
#      from this branch produces non-reproducible PDFs and is flagged.
DOC_DATE="$(awk '/^\*\*Date\.\*\*[[:space:]]+/ {
    s = $0
    sub(/^\*\*Date\.\*\*[[:space:]]+/, "", s)
    sub(/\.[[:space:]]*$/, "", s)
    print s; exit
}' "$SRC")"
if [ -z "$DOC_DATE" ]; then
    LATEST_TAG="$(git -C "$REPO_ROOT" tag --list 'memo-*' --sort=-creatordate 2>/dev/null | head -n1 || true)"
    if [ -n "$LATEST_TAG" ]; then
        DOC_DATE="$(echo "$LATEST_TAG" | sed -E 's/^memo-([0-9]{4}-[0-9]{2}-[0-9]{2}).*/\1/')"
    fi
fi
if [ -z "$DOC_DATE" ]; then
    DOC_DATE="$(date +%Y-%m-%d)"
    echo "WARNING: no **Date.** in $SRC and no memo-* tag on disk — falling back to wall-clock date $DOC_DATE; PDF will NOT be reproducible from this source." >&2
fi

# Export SOURCE_DATE_EPOCH from DOC_DATE so that Typst (and pandoc) embed a
# reproducible PDF CreationDate/ModDate. Without this, every re-render writes
# wall-clock timestamps into the PDF metadata, so two renders of the same
# tagged source produce PDFs that diff at the metadata level even when their
# rendered content is identical. With SOURCE_DATE_EPOCH set, Typst's
# embedded timestamps are deterministic and the resulting PDFs are
# byte-identical across re-renders. Verified: typst 0.14.2 + pandoc 3.9
# honor this. Date is interpreted at 00:00:00 UTC of DOC_DATE.
SOURCE_DATE_EPOCH="$(python3 -c "
import datetime, sys
y, m, d = (int(x) for x in '$DOC_DATE'.split('-'))
print(int(datetime.datetime(y, m, d, tzinfo=datetime.timezone.utc).timestamp()))
" 2>/dev/null)"
if [ -z "$SOURCE_DATE_EPOCH" ]; then
    echo "WARNING: could not derive SOURCE_DATE_EPOCH from DOC_DATE=$DOC_DATE — embedded PDF metadata will be wall-clock and PDF will not be byte-identical across re-renders." >&2
else
    export SOURCE_DATE_EPOCH
fi

# mktemp -t on macOS doesn't preserve the trailing .md, so pandoc can't
# deduce the format from extension. Create a directory and put a .md
# file inside it — portable and avoids the warning.
TMPDIR="$(mktemp -d -t pdf_render.XXXXXX)"
TMP="${TMPDIR}/src.md"
trap 'rm -rf "$TMPDIR"' EXIT

# Strip the leading title block from the H1 line up to (and including)
# the first horizontal rule `---` on its own line. Both PAPER.md and
# MANUSCRIPT.md use the convention of a `---` separator between the
# multi-paragraph title block and the first body section; the earlier
# "stop at the first blank line" awk only handled single-paragraph
# title blocks and leaked the second paragraph into the rendered body.
awk '
  BEGIN { in_meta = 0; stripped_h1 = 0 }
  !stripped_h1 && /^# / { stripped_h1 = 1; in_meta = 1; next }
  in_meta {
      if ($0 ~ /^---[[:space:]]*$/) { in_meta = 0 }
      next  # drop everything inside the title block, including the --- line itself
  }
  { print }
' "$SRC" > "$TMP"

# Resource path: figures referenced via relative paths
# (e.g. ![Figure 1](docs/figures/fig1...)) need to resolve from the
# repo root, not from /tmp where TMP lives. --resource-path tells
# pandoc where to look for those.
# `-f markdown-implicit_figures` stops pandoc from wrapping every solo
# image in a typst figure block. Otherwise typst auto-prefixes the
# alt text with "Figure N:", which collides with the markdown
# "**Figure N.**" bold caption paragraph immediately below — the
# rendered PDF would show "Figure N" three times per figure (typst
# auto-prefix + alt text + bold caption).
pandoc "$TMP" -o "$OUT" --pdf-engine=typst \
    -f markdown-implicit_figures \
    --resource-path="${REPO_ROOT}" \
    --metadata title="$TITLE" \
    --metadata subtitle="$SUBTITLE" \
    --metadata author="$AUTHOR" \
    --metadata date="$DOC_DATE" \
    --toc --toc-depth=2 \
    --variable=fontsize=10pt

ls -lh "$OUT"
