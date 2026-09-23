#!/usr/bin/env bash
set -euo pipefail

VALIDATION_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${VALIDATION_DIR}/../.." && pwd)"
export RUN_ROOT="${RUN_ROOT:-${1:-${ROOT_DIR}/validation_runs/reduced}}"

"${VALIDATION_DIR}/00_preflight.sh"
"${VALIDATION_DIR}/01_checkpoint_resume.sh"
"${VALIDATION_DIR}/02_cjejuni_matrix.sh"
"${VALIDATION_DIR}/03_compare.sh"
if [[ -n "${SECTION_D_REFERENCE:-}" ]]; then
    python3 "${VALIDATION_DIR}/check_section_d.py" \
        --run-root "$RUN_ROOT" --reference "$SECTION_D_REFERENCE" --image "$IMAGE"
fi

echo
echo 'Reduced same-source validation passed.'
