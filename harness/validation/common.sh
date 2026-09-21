#!/usr/bin/env bash
set -euo pipefail

VALIDATION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${VALIDATION_DIR}/../.." && pwd)"
IMAGE="${IMAGE:-nerve:v0.1.0}"
PSORT_CONTAINER="${PSORT_CONTAINER:-psortb}"
RUN_ROOT="${RUN_ROOT:-${ROOT_DIR}/validation_runs/reduced}"
INPUT_FASTA="${INPUT_FASTA:-${ROOT_DIR}/validation_data/campylobacter_jejuni/UP000000799_NCTC11168.fasta}"
CONSERVATION_FASTA="${CONSERVATION_FASTA:-${ROOT_DIR}/validation_data/campylobacter_jejuni/UP000000646_81-176.fasta}"

container_path() {
    case "$1" in
        "${ROOT_DIR}"/*) printf '/workdir/%s\n' "${1#"${ROOT_DIR}/"}" ;;
        *) printf '%s\n' "$1" ;;
    esac
}

require_new_dir() {
    if [ -e "$1" ]; then
        echo "Refusing to overwrite existing path: $1" >&2
        exit 1
    fi
    mkdir -p "$1"
}

run_nerve() {
    docker run --rm --network nerve-network \
        -v "${ROOT_DIR}:/workdir" \
        "${IMAGE}" "$@"
}

COMMON_ARGS=(
    --proteome1 "$(container_path "${INPUT_FASTA}")"
    --proteome2 "$(container_path "${CONSERVATION_FASTA}")"
    --gram n
    --annotation True
    --razor True
    --mouse True
    --virulent True
    --select True
    --epitopes True
    --mhci_length 9
    --mhcii_length 11
    --mhci_overlap 1
    --mhcii_overlap 1
    --epitope_percentile 0.9
    --ep_plots False
    --keep_snapshots True
)
