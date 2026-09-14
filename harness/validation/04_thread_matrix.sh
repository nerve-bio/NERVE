#!/usr/bin/env bash
set -euo pipefail
# Prove BLAST result stability at 1, an intermediate value, and available logical capacity.
: "${INPUT_FASTA:=validation_data/campylobacter_jejuni/UP000000799_NCTC11168.fasta}"
: "${CONSERVATION_FASTA:=validation_data/campylobacter_jejuni/UP000000646_81-176.fasta}"
for threads in 1 4 20; do
  ./NERVE.sh -p1 "$INPUT_FASTA" -p2 "$CONSERVATION_FASTA" -g n -t "$threads" -ew 1 \
    -wd "runs/threads-$threads" -ks True
done
for threads in 4 20; do
  python3 harness/compare_states.py --dir runs/threads-1/snapshots \
    "runs/threads-$threads/snapshots" --tol 0
  python3 harness/compare_candidates.py runs/threads-1/vaccine_candidates.csv \
    "runs/threads-$threads/vaccine_candidates.csv" --tol 0
done
