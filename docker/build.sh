#!/usr/bin/env bash
# Build the gx-nipt Docker image.
#
# Usage:
#   ./docker/build.sh                  # incremental (use layer cache)
#   ./docker/build.sh --no-cache       # full rebuild, ignore all cache
#   ./docker/build.sh --bump           # bump CACHE_BREAK and rebuild
#   ./docker/build.sh --no-cache --bump
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCKERFILE="${REPO_ROOT}/docker/Dockerfile"
IMAGE_NAME="gx-nipt:latest"

NO_CACHE=false
BUMP=false

for arg in "$@"; do
    case "$arg" in
        --no-cache) NO_CACHE=true ;;
        --bump)     BUMP=true ;;
        --help|-h)
            echo "Usage: $0 [--no-cache] [--bump]"
            echo "  --no-cache  Pass --no-cache to docker build (ignore all cached layers)"
            echo "  --bump      Increment CACHE_BREAK in Dockerfile before building"
            exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

# ── Bump CACHE_BREAK if requested ─────────────────────────────────────────────
if $BUMP; then
    current=$(grep -oP 'CACHE_BREAK=\K[0-9]+' "$DOCKERFILE")
    new=$((current + 1))
    sed -i "s/CACHE_BREAK=${current}/CACHE_BREAK=${new}/" "$DOCKERFILE"
    echo "[build.sh] CACHE_BREAK: ${current} → ${new}"
fi

# ── Build ──────────────────────────────────────────────────────────────────────
BUILD_ARGS=()
$NO_CACHE && BUILD_ARGS+=(--no-cache)

echo "[build.sh] Building ${IMAGE_NAME} ..."
echo "[build.sh] Context : ${REPO_ROOT}"
echo "[build.sh] Options : ${BUILD_ARGS[*]:-none}"
echo ""

START=$(date +%s)

docker build \
    "${BUILD_ARGS[@]}" \
    -t "${IMAGE_NAME}" \
    -f "${DOCKERFILE}" \
    "${REPO_ROOT}"

END=$(date +%s)
ELAPSED=$(( END - START ))
echo ""
echo "[build.sh] Done in ${ELAPSED}s → ${IMAGE_NAME}"
docker images "${IMAGE_NAME%:*}" --format "  ID={{.ID}}  Size={{.Size}}  Created={{.CreatedAt}}"
