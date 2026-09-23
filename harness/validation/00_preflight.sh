#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

cd "$ROOT_DIR"
git diff --check
python3 -m compileall -q code harness tests
bash -n NERVE.sh harness/validation/*.sh

test -s "$INPUT_FASTA" && test -s "$CONSERVATION_FASTA" || {
    echo 'The fixed C. jejuni input FASTAs are missing or empty.' >&2
    exit 1
}
test "$(sha256sum "$INPUT_FASTA" | cut -d' ' -f1)" = \
    'baf55425a1477aa534de8c21868cc2fba8bee8f66d156234a7b9233708394df7' || {
    echo "Unexpected C. jejuni primary input: $INPUT_FASTA" >&2
    exit 1
}
test "$(sha256sum "$CONSERVATION_FASTA" | cut -d' ' -f1)" = \
    'fbed332360663d26c34b478c96884370dfc38e61696299adf5033137670bed18' || {
    echo "Unexpected C. jejuni comparison input: $CONSERVATION_FASTA" >&2
    exit 1
}

IMAGE_VERSION="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.version" }}' "$IMAGE")"
IMAGE_TREE="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.source-tree" }}' "$IMAGE")"
SOURCE_TREE="$(git ls-files -co --exclude-standard code Dockerfile .dockerignore NERVE.sh | \
    grep -Ev '(^|/)(__pycache__/|.*\.pyc$)' | LC_ALL=C sort | \
    xargs sha256sum | sha256sum | cut -d' ' -f1)"

test "$IMAGE_VERSION" = "v0.1.0"
test "$IMAGE_TREE" = "$SOURCE_TREE"
test "$(docker inspect -f '{{.State.Running}}' "$PSORT_CONTAINER" 2>/dev/null)" = "true" || {
    echo "PSORTb container \"${PSORT_CONTAINER}\" must be running." >&2
    exit 1
}
test "$(docker inspect -f '{{if index .NetworkSettings.Networks "nerve-network"}}true{{end}}' "$PSORT_CONTAINER")" = "true" || {
    echo "PSORTb container \"${PSORT_CONTAINER}\" must use nerve-network." >&2
    exit 1
}

docker run --rm \
    -v "${ROOT_DIR}/tests:/usr/nerve_python/NERVE/tests:ro" \
    -v "${ROOT_DIR}/harness:/usr/nerve_python/NERVE/harness:ro" \
    -w /usr/nerve_python/NERVE \
    --entrypoint python3 "$IMAGE" \
    -m pytest -q tests/test_equivalence.py tests/test_harness.py -p no:debugging

echo "Preflight passed for ${IMAGE}."
