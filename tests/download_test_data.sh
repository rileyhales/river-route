#!/usr/bin/env bash
set -euxo pipefail
# Downloads test data for the river-route test suite.

# the zip contains the directories:
# era5 - containing 12 example months of era5 data
# discharge
# |- vpu=208
# |  |- discharge_*.nc
# |- vpu=410
# |  |- discharge_*.nc
# |- vpu=718
#    |- discharge_*.nc
# qlateral
# |- vpu=208
# |  |- qlateral_*.nc
# |- vpu=410
# |  |- qlateral_*.nc
# |- vpu=718
#    |- qlateral_*.nc

# first download era5 and correct qlateral and discharges
wget -q --tries=3 -O routing-test-data.zip "https://geoglows-v2.s3.amazonaws.com/routing-test-data.zip"
# then unzip it to the test data folder and remove the zip
unzip -qo routing-test-data.zip -d tests/data
rm -f routing-test-data.zip

# first see what the list of vpu=* directories are contained in tests/data/discharge
VPUS=$(find tests/data/discharge -type d -name "vpu=*")

# for each vpu that test data is provided for, download the hydrography and routing config inputs
for VPU in $VPUS; do
    VPU_NAME=$(basename "${VPU}")
    echo "Downloading test data for ${VPU_NAME} ..."
    s5cmd --no-sign-request cp "s3://geoglows-v2/routing-configs/${VPU_NAME}/*" tests/data/routing-configs/${VPU_NAME}/
    s5cmd --no-sign-request cp "s3://geoglows-v2/hydrography/${VPU_NAME}/*" tests/data/hydrography/${VPU_NAME}/
done

# the inputs then need to be processed by the test suite which happens in another script
