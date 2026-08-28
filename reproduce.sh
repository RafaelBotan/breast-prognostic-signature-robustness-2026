#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${COHORT_DIR:-}" ]]; then
  echo "Set COHORT_DIR to the local processed public-cohort directory." >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
export STRESS_PROJECT_DIR="$REPO_ROOT"
export SUBMITTED_AGGREGATE_DIR="$REPO_ROOT/data/reference_aggregate"
export REVISION_OUTPUT_DIR="${REVISION_OUTPUT_DIR:-$REPO_ROOT/outputs}"
export REVISION_CACHE_DIR="${REVISION_CACHE_DIR:-$REPO_ROOT/cache}"
export FIGURE_OUTPUT_DIR="${FIGURE_OUTPUT_DIR:-$REPO_ROOT/outputs/figures}"
export B_MATCH="${B_MATCH:-1000}"

mkdir -p "$REVISION_OUTPUT_DIR" "$REVISION_CACHE_DIR" "$FIGURE_OUTPUT_DIR"

Rscript "$REPO_ROOT/code/revision/reviewer1_revision_analysis.R"
Rscript "$REPO_ROOT/code/revision/postprocess_matched_outputs.R"
Rscript "$REPO_ROOT/code/revision/stress_e2_platform_revision.R"
Rscript "$REPO_ROOT/code/revision/make_figure1_revision.R"
Rscript "$REPO_ROOT/code/revision/make_figure2_revision.R"

echo "Revision outputs written to $REVISION_OUTPUT_DIR"

