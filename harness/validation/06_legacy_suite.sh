#!/usr/bin/env bash
set -euo pipefail
# Last gate: run the repository's older Docker-oriented cases against the verified port image.
docker run --rm -v "$PWD:/workdir" --entrypoint python3 nerve:v0.0.8 \
  -m pytest -q /workdir/tests/tests.py
# TODO convert implicit/manual success criteria in tests.py into assertions and nonzero failures.
