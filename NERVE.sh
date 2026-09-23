#!/usr/bin/env bash
set -euo pipefail

PSORT_IMAGE='francecosta/psortb_http_api@sha256:68357a19efee3aa837f94e3af84c19ad8e77bf19b86d64a9a22282a5447d6ea8'
BASE_IMAGE='francecosta/nerve@sha256:f49feae627370f12bde8a4ec12a117b53e5840caf10a21531e7fc1e9c46fd9f9'
IMAGE='nerve:v0.1.0'
USE_IMAGE=false
BUILD_ONLY=false
INTERACTIVE=false
while (( $# )); do
    case "$1" in
        --build-only) BUILD_ONLY=true; shift ;;
        --image)
            if (( $# < 2 )); then
                echo '--image requires an image ID or digest' >&2
                exit 2
            fi
            IMAGE="$2"; USE_IMAGE=true; shift 2 ;;
        --interactive) INTERACTIVE=true; shift ;;
        --) shift; break ;;
        *) break ;;
    esac
done

SOURCE_COMMIT="$(git rev-parse HEAD)"
SOURCE_TREE="$(git ls-files -co --exclude-standard code Dockerfile .dockerignore NERVE.sh | \
    grep -Ev '(^|/)(__pycache__/|.*\.pyc$)' | LC_ALL=C sort | \
    xargs sha256sum | sha256sum | cut -d' ' -f1)"
if git diff --quiet --ignore-submodules HEAD -- code Dockerfile .dockerignore NERVE.sh && \
   [ -z "$(git ls-files --others --exclude-standard code Dockerfile .dockerignore NERVE.sh)" ]; then
    SOURCE_DIRTY=false
else
    SOURCE_DIRTY=true
fi

if $BUILD_ONLY && $USE_IMAGE; then
    echo '--build-only cannot be combined with --image' >&2
    exit 2
fi
if ! $USE_IMAGE; then
    echo "Building $IMAGE from $SOURCE_COMMIT (dirty=$SOURCE_DIRTY)"
    docker build --build-arg NERVE_SOURCE_COMMIT="$SOURCE_COMMIT" \
        --build-arg NERVE_SOURCE_DIRTY="$SOURCE_DIRTY" \
        --build-arg NERVE_SOURCE_TREE="$SOURCE_TREE" -t "$IMAGE" .
fi

IMAGE_COMMIT="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.revision" }}' "$IMAGE")"
IMAGE_TREE="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.source-tree" }}' "$IMAGE")"
if [[ "$IMAGE_COMMIT" != "$SOURCE_COMMIT" || "$IMAGE_TREE" != "$SOURCE_TREE" ]]; then
    echo "Image source mismatch: image=$IMAGE_COMMIT/$IMAGE_TREE checkout=$SOURCE_COMMIT/$SOURCE_TREE" >&2
    exit 1
fi
IMAGE_ID="$(docker image inspect -f '{{.Id}}' "$IMAGE")"
echo "Verified image source: commit=$IMAGE_COMMIT tree=$IMAGE_TREE image=$IMAGE_ID"
$BUILD_ONLY && exit 0

if ! docker network inspect nerve-network >/dev/null 2>&1; then
    docker network create nerve-network >/dev/null
fi
if docker container inspect psortb >/dev/null 2>&1; then
    if [[ "$(docker inspect -f '{{.State.Running}}' psortb)" != true ]] || \
       [[ "$(docker inspect -f '{{if index .NetworkSettings.Networks "nerve-network"}}true{{end}}' psortb)" != true ]]; then
        echo 'Existing psortb container must be running on nerve-network.' >&2
        exit 1
    fi
else
    docker run --rm -d --name psortb -p 8080:8080 --network nerve-network "$PSORT_IMAGE" >/dev/null
fi
PSORT_ID="$(docker inspect -f '{{.Image}}' psortb)"

TTY=()
if $INTERACTIVE; then
    if [[ ! -t 0 ]]; then
        echo '--interactive requires a terminal on stdin' >&2
        exit 2
    fi
    TTY=(-it)
fi
docker run --rm --network nerve-network -p 8880:8880 "${TTY[@]}" \
    -e NERVE_IMAGE_ID="$IMAGE_ID" \
    -e NERVE_BASE_IMAGE_DIGEST="$BASE_IMAGE" \
    -e NERVE_PSORTB_IMAGE_ID="$PSORT_ID" \
    -v "$(pwd):/workdir" "$IMAGE" "$@"
