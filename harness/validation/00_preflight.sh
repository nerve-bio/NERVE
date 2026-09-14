#!/usr/bin/env bash
set -euo pipefail
# Fast gate: source/image identity, syntax, and focused regression tests.
./NERVE.sh --build-only
docker run --rm -v "$PWD:/workdir" --entrypoint python3 nerve:v0.0.8 \
  -m pytest -q /workdir/tests/test_equivalence.py
docker image inspect nerve:v0.0.8
