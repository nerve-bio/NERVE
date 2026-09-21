#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

MATRIX_ROOT="${RUN_ROOT}/cjejuni"
REFERENCE_DIR="${MATRIX_ROOT}/reference-t1-ew1"
BLAST_DIR="${MATRIX_ROOT}/blast-t20"
EPITOPE_DIR="${MATRIX_ROOT}/epitope-ew20"

compare_state() {
    python3 "${ROOT_DIR}/harness/compare_states.py" "$1" "$2" \
        --attrs "$3" --tol 0
}

compare_state "${REFERENCE_DIR}/snapshots/06_autoimmunity.pkl" \
    "${BLAST_DIR}/06_autoimmunity.pkl" \
    'accession,sequence,list_of_shared_human_peps,sapiens_peptides_sum,list_of_peptides_from_comparison_with_mhcpep_sapiens'
compare_state "${REFERENCE_DIR}/snapshots/07_mouse.pkl" \
    "${BLAST_DIR}/07_mouse.pkl" \
    'accession,sequence,list_of_shared_mouse_peps,mouse_peptides_sum,list_of_peptides_from_comparison_with_mhcpep_mouse'
compare_state "${REFERENCE_DIR}/snapshots/08_conservation.pkl" \
    "${BLAST_DIR}/08_conservation.pkl" \
    'accession,sequence,list_of_shared_conserv_proteome_peps,conservation_score'

compare_state "${REFERENCE_DIR}/snapshots/12_epitope.pkl" \
    "${EPITOPE_DIR}/12_epitope.pkl" \
    'accession,score,MHC1_binders,MHC2_binders,MHC1_pb_binders,MHC2_pb_binders'

python3 "${ROOT_DIR}/harness/validation/compare_epitope_outputs.py" \
    "${REFERENCE_DIR}/epitope" "${EPITOPE_DIR}/epitope" --tol 0

echo 'C. jejuni BLAST-thread and epitope-worker comparisons passed exactly.'
