#! /bin/bash -e
# script to run NERVE standalone trough docker

PSORT_VERSION="v0.0.1"
NERVE_VERSION="v0.1.0"

# Record the source revision in the image.
NERVE_SOURCE_COMMIT="$(git rev-parse HEAD 2>/dev/null || printf unknown)"
NERVE_SOURCE_TREE="$(git ls-files -co --exclude-standard code Dockerfile .dockerignore NERVE.sh | \
    grep -Ev '(^|/)(__pycache__/|.*\.pyc$)' | LC_ALL=C sort | \
    xargs sha256sum | sha256sum | cut -d' ' -f1)"
if git diff --quiet --ignore-submodules HEAD -- 2>/dev/null && \
   [ -z "$(git ls-files --others --exclude-standard 2>/dev/null)" ]; then
    NERVE_SOURCE_DIRTY="false"
else
    NERVE_SOURCE_DIRTY="true"
fi

# Build the current source.
echo "Building nerve:${NERVE_VERSION} from ${NERVE_SOURCE_COMMIT} (dirty=${NERVE_SOURCE_DIRTY})"
docker build \
    --build-arg NERVE_SOURCE_COMMIT="${NERVE_SOURCE_COMMIT}" \
    --build-arg NERVE_SOURCE_DIRTY="${NERVE_SOURCE_DIRTY}" \
    --build-arg NERVE_SOURCE_TREE="${NERVE_SOURCE_TREE}" \
    -t "nerve:${NERVE_VERSION}" .

IMAGE_COMMIT="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.revision" }}' "nerve:${NERVE_VERSION}")"
IMAGE_TREE="$(docker image inspect -f '{{ index .Config.Labels "org.opencontainers.image.source-tree" }}' "nerve:${NERVE_VERSION}")"
if [ "${IMAGE_COMMIT}" != "${NERVE_SOURCE_COMMIT}" ] || [ "${IMAGE_TREE}" != "${NERVE_SOURCE_TREE}" ]; then
    echo "Image provenance verification failed" >&2
    exit 1
fi
echo "Verified image source: commit=${IMAGE_COMMIT} tree=${IMAGE_TREE}"

if [ "${1:-}" = "--build-only" ]; then
    exit 0
fi

# Start the runtime services after the build.
[ ! "$(docker network ls | grep nerve-network)" ] && docker network create nerve-network --attachable
[ ! "$(docker ps | grep francecosta/psortb_http_api:$PSORT_VERSION)" ] && \
    docker run --rm -p 8080:8080 --network nerve-network --name psortb -d \
    "francecosta/psortb_http_api:$PSORT_VERSION"

# run nerve container
docker run --network nerve-network -p 8880:8880 -it -v "$(pwd):/workdir" "nerve:${NERVE_VERSION}" "$@"
