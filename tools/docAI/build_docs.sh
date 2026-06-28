#!/usr/bin/env bash
# Build both the standard and AI-enhanced Sphinx documentation.
# Run from the repository root.
#
# Usage:
#   ./tools/docAI/build_docs.sh [options] [OUTPUT_DIR]
#
# Options:
#   --skip-docai   Skip the doc_ai.py LLM step (use existing ai_content/ tree)
#   --model MODEL  LLM model for doc_ai.py (default: gpt-5.4-nano)
#
# OUTPUT_DIR defaults to tools/docAI/_build/html.
# Produces:
#   OUTPUT_DIR/standard/   — standard documentation
#   OUTPUT_DIR/ai/         — AI-enhanced documentation

set -euo pipefail

SKIP_DOCAI=0
MODEL="gpt-5.4-nano"
OUTPUT_DIR="tools/docAI/_build/html"
AI_CONTENT_DIR="tools/docAI/ai_content"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-docai) SKIP_DOCAI=1; shift ;;
        --model)      MODEL="$2"; shift 2 ;;
        *)            OUTPUT_DIR="$1"; shift ;;
    esac
done

echo "=== Step 1: generate C++ API RST (standard) ==="
python3 tools/docAI/generate_cpp_docs.py

echo "=== Step 2: generate tutorials RST ==="
python3 tools/docAI/generate_tutorials_docs.py

echo "=== Step 3: Sphinx standard build ==="
# -E forces a fresh environment each run. Needed because the AI build injects
# ai_content via a source-read hook that Sphinx does not track for cache
# invalidation: an incremental build would silently miss ai_content-only edits.
# Kept on the standard build too so both variants are always consistent.
sphinx-build -E -b html -j auto documentation/sphinx "${OUTPUT_DIR}/standard"

if [[ $SKIP_DOCAI -eq 0 ]]; then
    echo "=== Step 4: run doc_ai to populate AI content tree ==="
    python3 tools/docAI/doc_ai.py \
        --mode all \
        --model "${MODEL}" \
        --output-dir "${AI_CONTENT_DIR}"
else
    echo "=== Step 4: skipped doc_ai (--skip-docai) ==="
fi

echo "=== Step 5: generate C++ API RST (AI-enhanced) ==="
python3 tools/docAI/generate_cpp_docs.py --ai-content-dir "${AI_CONTENT_DIR}"

echo "=== Step 6: Sphinx AI-enhanced build ==="
# -E (fresh environment) is required here: the source-read overlay that injects
# ai_content is not cache-tracked, so without it an incremental rebuild keeps
# stale doctrees and the AI edits do not appear.
SPHINX_AI=1 sphinx-build -E -b html -j auto documentation/sphinx "${OUTPUT_DIR}/ai"

echo ""
echo "Done. Output:"
echo "  Standard : ${OUTPUT_DIR}/standard"
echo "  AI        : ${OUTPUT_DIR}/ai"
