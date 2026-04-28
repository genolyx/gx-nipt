#!/usr/bin/env bash
# =========================================================
#  gx-nipt build helper
#
#  Usage:
#    ./build.sh                          # build gx-nipt:latest (cached)
#    ./build.sh --no-cache               # full rebuild, ignore all cache
#    ./build.sh --bump                   # bump CACHE_BREAK and rebuild
#    ./build.sh --tag 1.2.0              # also tag as gx-nipt:1.2.0
#    ./build.sh --tag 1.2.0 --push       # tag + push to genolyx/ registry
#    ./build.sh --push                   # push :latest only
# =========================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REGISTRY="genolyx"
IMAGE="gx-nipt"
LOCAL_TAG="latest"

VERSION=""
PUSH=false
PASSTHROUGH_ARGS=()

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag|-t)       VERSION="$2"; shift 2 ;;
        --push)         PUSH=true; shift ;;
        --no-cache)     PASSTHROUGH_ARGS+=(--no-cache); shift ;;
        --bump)         PASSTHROUGH_ARGS+=(--bump); shift ;;
        --help|-h)
            sed -n '3,10p' "${BASH_SOURCE[0]}" | sed 's/^#  //'
            exit 0 ;;
        *) echo "[build.sh] Unknown option: $1" >&2; exit 1 ;;
    esac
done

LOCAL_IMAGE="${IMAGE}:${LOCAL_TAG}"
REGISTRY_LATEST="${REGISTRY}/${IMAGE}:latest"

# ── 1. Build local image ──────────────────────────────────────────────────────
"${REPO_ROOT}/docker/build.sh" "${PASSTHROUGH_ARGS[@]+"${PASSTHROUGH_ARGS[@]}"}"

# ── 2. Version tag (local) ───────────────────────────────────────────────────
if [[ -n "$VERSION" ]]; then
    docker tag "${LOCAL_IMAGE}" "${IMAGE}:${VERSION}"
    echo "[build.sh] Tagged  : ${IMAGE}:${VERSION}"
fi

# ── 3. Push to registry ───────────────────────────────────────────────────────
if $PUSH; then
    docker tag "${LOCAL_IMAGE}" "${REGISTRY_LATEST}"
    echo "[build.sh] Pushing : ${REGISTRY_LATEST}"
    docker push "${REGISTRY_LATEST}"

    if [[ -n "$VERSION" ]]; then
        REGISTRY_VERSION="${REGISTRY}/${IMAGE}:${VERSION}"
        docker tag "${LOCAL_IMAGE}" "${REGISTRY_VERSION}"
        echo "[build.sh] Pushing : ${REGISTRY_VERSION}"
        docker push "${REGISTRY_VERSION}"
    fi

    echo "[build.sh] Push complete"
fi
