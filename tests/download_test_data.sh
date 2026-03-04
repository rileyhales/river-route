#!/usr/bin/env bash
# Downloads test data for the river-route test suite.
# Data is hosted as a compressed archive on cloud storage.
# Set DATA_URL environment variable to override the default download location.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
ARCHIVE="${SCRIPT_DIR}/.test-data.zip"

# ── Configuration ────────────────────────────────────────────────────────────
DEFAULT_URL="https://geoglows-v2.s3.amazonaws.com/routing-test-data.zip"
DATA_URL="${DATA_URL:-$DEFAULT_URL}"

# ── Skip if data already exists ──────────────────────────────────────────────
if [ -f "${DATA_DIR}/routing-configs/vpu=718/params.parquet" ] && \
   [ -f "${DATA_DIR}/era5/era5_194001.nc" ]; then
    echo "Test data already present, skipping download."
    exit 0
fi

# ── Download ─────────────────────────────────────────────────────────────────
echo "Downloading test data from ${DATA_URL} ..."
if command -v curl &>/dev/null; then
    curl -fSL --retry 3 --retry-delay 5 -o "${ARCHIVE}" "${DATA_URL}"
elif command -v wget &>/dev/null; then
    wget -q --tries=3 -O "${ARCHIVE}" "${DATA_URL}"
else
    echo "Error: neither curl nor wget found" >&2
    exit 1
fi

# ── Extract ──────────────────────────────────────────────────────────────────
echo "Extracting test data to ${DATA_DIR} ..."
mkdir -p "${DATA_DIR}"
unzip -qo "${ARCHIVE}" -d "${DATA_DIR}"
rm -f "${ARCHIVE}"

echo "Test data ready."
