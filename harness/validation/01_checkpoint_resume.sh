#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

CHECK_ROOT="${RUN_ROOT}/checkpoint"
RESUME_DIR="${CHECK_ROOT}/resume"
MISMATCH_DIR="${CHECK_ROOT}/mismatch"
FRESH_DIR="${CHECK_ROOT}/fresh"
require_new_dir "$CHECK_ROOT"
mkdir -p "$RESUME_DIR" "$MISMATCH_DIR" "$FRESH_DIR"

CHECK_INPUT="${CHECK_ROOT}/input_100.fasta"
awk '/^>/{n++} n>100{exit} {print}' "$INPUT_FASTA" > "$CHECK_INPUT"
test "$(grep -c '^>' "$CHECK_INPUT")" -eq 100
SMALL_FASTA="$(container_path "$CHECK_INPUT")"
CHECK_ARGS=(
    --proteome1 "$SMALL_FASTA" --gram n
    --annotation False --razor False --mouse False --virulent False
    --select True --epitopes False --ep_plots False
    --threads 1 --epitope_workers 1 --keep_snapshots True
)

python3 "${VALIDATION_DIR}/run_checkpoint_interrupt.py" \
    --root "$ROOT_DIR" --image "$IMAGE" --output "$RESUME_DIR" --input "$CHECK_INPUT"

test -s "${RESUME_DIR}/checkpoint.pkl"
test ! -e "${RESUME_DIR}/checkpoint.pkl.tmp"
cp "${RESUME_DIR}/checkpoint.pkl" "${MISMATCH_DIR}/checkpoint.pkl"
cp "${RESUME_DIR}/run_manifest.json" "${MISMATCH_DIR}/run_manifest.json"

MISMATCH_LOG="${MISMATCH_DIR}/expected-rejection.log"
set +e
run_nerve "${CHECK_ARGS[@]}" --gram p --resume True \
    --working_dir "$(container_path "$MISMATCH_DIR")" \
    >"$MISMATCH_LOG" 2>&1
MISMATCH_STATUS=$?
set -e
test "$MISMATCH_STATUS" -ne 0
grep -q 'Checkpoint does not match' "$MISMATCH_LOG"
echo 'Incompatible checkpoint rejected.'

run_nerve "${CHECK_ARGS[@]}" --resume True \
    --working_dir "$(container_path "$RESUME_DIR")"

test -s "${RESUME_DIR}/vaccine_candidates.csv"
test ! -e "${RESUME_DIR}/checkpoint.pkl"
test -f "${RESUME_DIR}/snapshots/11_select.pkl"

run_nerve "${CHECK_ARGS[@]}" \
    --working_dir "$(container_path "$FRESH_DIR")"
python3 "${ROOT_DIR}/harness/compare_states.py" --dir \
    "${FRESH_DIR}/snapshots" "${RESUME_DIR}/snapshots" --tol 0 \
    --canonical-localization-ties
python3 "${ROOT_DIR}/harness/compare_candidates.py" \
    "${FRESH_DIR}/vaccine_candidates.csv" \
    "${RESUME_DIR}/vaccine_candidates.csv" --tol 0
echo 'Checkpoint interruption, mismatch rejection and exact fresh-run comparison passed.'
