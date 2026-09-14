#!/usr/bin/env bash
set -euo pipefail
# Repeat baseline/port comparison on adversarial headers, but whitelist only intended ID fixes.
: "${REFERENCE_DIR:?set completed adversarial-header baseline directory}"
: "${PORT_DIR:?set completed adversarial-header port directory}"
python3 harness/compare_states.py --dir "$REFERENCE_DIR/snapshots" "$PORT_DIR/snapshots" --tol 0 \
  > "$PORT_DIR/id_delta_report.txt" || true
# TODO machine-readable allowlist: canonical-header changes, corrected exact protein assignment,
# and their downstream consequences. Fail for every delta not reachable from that allowlist.
