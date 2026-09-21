#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

MATRIX_ROOT="${RUN_ROOT}/cjejuni"
REFERENCE_DIR="${MATRIX_ROOT}/reference-t1-ew1"
BLAST_DIR="${MATRIX_ROOT}/blast-t20"
EPITOPE_DIR="${MATRIX_ROOT}/epitope-ew20"
require_new_dir "$MATRIX_ROOT"

run_nerve "${COMMON_ARGS[@]}" \
    --threads 1 --epitope_workers 1 \
    --working_dir "$(container_path "$REFERENCE_DIR")"

docker run --rm --network nerve-network \
    -v "${ROOT_DIR}:/workdir" --entrypoint python3 "$IMAGE" \
    /workdir/harness/validation/run_blast_variant.py \
    --snapshot "$(container_path "$REFERENCE_DIR")/snapshots/05_razor.pkl" \
    --proteome1 "$(container_path "$REFERENCE_DIR")/cleaned_$(basename "$INPUT_FASTA")" \
    --proteome2 "$(container_path "$REFERENCE_DIR")/cleaned_$(basename "$CONSERVATION_FASTA")" \
    --output "$(container_path "$BLAST_DIR")" --threads 20

docker run --rm \
    -v "${ROOT_DIR}:/workdir" --entrypoint python3 "$IMAGE" \
    /workdir/harness/validation/run_epitope_variant.py \
    --snapshot "$(container_path "$REFERENCE_DIR")/snapshots/11_select.pkl" \
    --output "$(container_path "$EPITOPE_DIR")" --workers 20

echo "C. jejuni reference and phase variants completed in ${MATRIX_ROOT}."
