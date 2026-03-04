#!/usr/bin/env bash
# Packages tests/data/ into a tarball suitable for uploading to cloud storage.
# Run from the repo root:  ./tests/package_test_data.sh
#
# After running, upload the resulting archive to your bucket and set the
# DATA_URL (or GitHub Actions variable TEST_DATA_URL) to the public URL.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT="${SCRIPT_DIR}/routing-test-data.zip"

if [ ! -d "${DATA_DIR}" ]; then
    echo "Error: ${DATA_DIR} does not exist" >&2
    exit 1
fi

echo "Packaging test data ..."
cd "${DATA_DIR}"
zip -rq "${OUTPUT}" \
    era5 \
    routing-configs \
    discharge \
    catchment-volumes \
    hydrography \
    -x '*.DS_Store' '*__pycache__*' 'scripts/*'

SIZE=$(du -h "${OUTPUT}" | cut -f1)
echo "Created ${OUTPUT} (${SIZE})"
echo ""
echo "Next steps:"
echo "  1. Upload to your cloud bucket (S3, GCS, etc.)"
echo "  2. Set the public URL as GitHub Actions variable TEST_DATA_URL"
echo "     (Settings > Secrets and variables > Actions > Variables)"
echo "  3. Or set DATA_URL env var when running tests/download_test_data.sh locally"
