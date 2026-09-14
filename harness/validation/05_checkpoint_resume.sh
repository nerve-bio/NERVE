#!/usr/bin/env bash
set -euo pipefail
# Interrupt after at least one checkpoint, resume it, then compare with an uninterrupted run.
: "${INPUT_FASTA:=validation_data/campylobacter_jejuni/UP000000799_NCTC11168.fasta}"
: "${CONSERVATION_FASTA:=validation_data/campylobacter_jejuni/UP000000646_81-176.fasta}"
./NERVE.sh -p1 "$INPUT_FASTA" -p2 "$CONSERVATION_FASTA" -g n -wd runs/uninterrupted -ks True
timeout --signal=TERM 120 ./NERVE.sh -p1 "$INPUT_FASTA" -p2 "$CONSERVATION_FASTA" -g n \
  -wd runs/resumed -ks True -ck True || test $? -eq 124
test -s runs/resumed/checkpoint.pkl
./NERVE.sh -p1 "$INPUT_FASTA" -p2 "$CONSERVATION_FASTA" -g n -wd runs/resumed -ks True -rs True
python3 harness/compare_states.py --dir runs/uninterrupted/snapshots runs/resumed/snapshots --tol 0
python3 harness/compare_candidates.py runs/uninterrupted/vaccine_candidates.csv \
  runs/resumed/vaccine_candidates.csv --tol 0
test ! -e runs/resumed/checkpoint.pkl
