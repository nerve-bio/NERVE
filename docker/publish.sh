#! /bin/bash -e
# Build, push and record the image downstream consumers pin.
#   docker/publish.sh             build, test, push, rewrite the published block
#   docker/publish.sh --dry-run   build and test only
#   docker/publish.sh --archive   also keep a local archive of the pushed image

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOCK="$ROOT/base-image.lock"
REPO="ghcr.io/nerve-bio/nerve"
ARCHIVE_DIR="${NERVE_BASE_ARCHIVE_DIR:-$HOME/.cache/nerve-base}"

dry=0
archive=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) dry=1 ;;
        --archive) archive=1 ;;
        *) echo "usage: $0 [--dry-run] [--archive]" >&2; exit 1 ;;
    esac
done

# The published image names a commit, so the tree it is built from must be that commit.
dirty="$(git -C "$ROOT" status --porcelain=v1 -- code Dockerfile .dockerignore NERVE.sh tests)"
[ -z "$dirty" ] || { echo "tracked build or test source is dirty:" >&2; echo "$dirty" >&2; exit 1; }
COMMIT="$(git -C "$ROOT" rev-parse --short HEAD)"
IMAGE="$REPO:$COMMIT"

"$ROOT/docker/verify_base.sh"
docker build --build-arg NERVE_COMMIT="$COMMIT" -t "$IMAGE" "$ROOT"

# deploy.sh reads this label back off the pinned image
revision="$(docker image inspect "$IMAGE" --format '{{index .Config.Labels "org.opencontainers.image.revision"}}')"
[ "$revision" = "$COMMIT" ] || { echo "image revision $revision, expected $COMMIT" >&2; exit 1; }
docker run --rm -v "$ROOT":/repo -w /repo --entrypoint pytest "$IMAGE" -q tests/test_checks.py

[ "$dry" = 0 ] || { echo "dry run: $IMAGE built and tested, nothing pushed"; exit 0; }
docker push "$IMAGE"

CONFIG="$(docker image inspect "$IMAGE" --format '{{.Id}}')"
DIGEST="$(docker image inspect "$IMAGE" --format '{{range .RepoDigests}}{{println .}}{{end}}' \
    | grep -o 'sha256:[0-9a-f]*' | head -1)"
[ -n "$DIGEST" ] || { echo "no registry digest for $IMAGE; the push did not land" >&2; exit 1; }

ARCHIVE_NAME=""
ARCHIVE_SHA=""
if [ "$archive" = 1 ]; then
    mkdir -p "$ARCHIVE_DIR"
    ARCHIVE_NAME="nerve-bio-$COMMIT.tar"
    docker save -o "$ARCHIVE_DIR/$ARCHIVE_NAME" "$IMAGE"
    ARCHIVE_SHA="$(sha256sum "$ARCHIVE_DIR/$ARCHIVE_NAME" | cut -d' ' -f1)"
    echo "$ARCHIVE_SHA  $ARCHIVE_NAME" >> "$ARCHIVE_DIR/sha256.txt"
fi

# Absent archive fields are read as empty; verify_base.sh then restores by digest.
awk -v image="$IMAGE" -v digest="$DIGEST" -v config="$CONFIG" -v commit="$COMMIT" \
    -v name="$ARCHIVE_NAME" -v asha="$ARCHIVE_SHA" '
    $1 == "published_image"           { printf "%-22s %s\n", $1, image;  seen=1; next }
    $1 == "published_digest"          { printf "%-22s %s\n", $1, digest; next }
    $1 == "published_config"          { printf "%-22s %s\n", $1, config; next }
    $1 == "published_code_commit"     { printf "%-22s %s\n", $1, commit; next }
    $1 == "published_archive_file"    { if (name) printf "%-25s %s\n", $1, name; next }
    $1 == "published_archive_sha256"  { if (asha) printf "%-25s %s\n", $1, asha; next }
    { print }
    END { if (!seen) { print "no published_image line in the lock" > "/dev/stderr"; exit 1 } }
' "$LOCK" > "$LOCK.new"
mv "$LOCK.new" "$LOCK"

"$ROOT/docker/verify_base.sh"
echo
echo "pushed   $IMAGE"
echo "digest   $DIGEST"
echo "pin this in the consumer's worker base:"
echo "  FROM $REPO@$DIGEST AS base"
