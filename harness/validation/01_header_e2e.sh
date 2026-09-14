#!/usr/bin/env bash
set -euo pipefail
# Real pipeline fixture: slash-bearing, bare, duplicate, RefSeq and UniProt headers.
: "${HEADER_FASTA:?set HEADER_FASTA to the mixed-header fixture}"
: "${RUN_ROOT:=runs/headers}"
./NERVE.sh -p1 "$HEADER_FASTA" -g n -wd "$RUN_ROOT" -ks True
# TODO assertions: canonical accessions unique; no nested epitope paths; original ids retained;
# every kept record reaches external-tool lookup; QC audit rows explain every modification/discard.
