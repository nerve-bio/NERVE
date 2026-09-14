#!/usr/bin/env bash
set -euo pipefail
# One ordinary-header proteome, identical effective parameters, full reference and port runs.
: "${INPUT_FASTA:=validation_data/campylobacter_jejuni/UP000000799_NCTC11168.fasta}"
: "${CONSERVATION_FASTA:=validation_data/campylobacter_jejuni/UP000000646_81-176.fasta}"
: "${REFERENCE_RUNNER:?set command that launches the pinned baseline image}"
: "${PORT_RUNNER:=./NERVE.sh}"
common=(-p1 "$INPUT_FASTA" -p2 "$CONSERVATION_FASTA" -g n -ks True)
eval "$REFERENCE_RUNNER" "${common[@]}" -wd runs/reference
eval "$PORT_RUNNER" "${common[@]}" -wd runs/port
python3 harness/compare_states.py --dir runs/reference/snapshots runs/port/snapshots --tol 0
python3 harness/compare_candidates.py runs/reference/vaccine_candidates.csv runs/port/vaccine_candidates.csv --tol 0
diff -ru runs/reference/epitope runs/port/epitope
