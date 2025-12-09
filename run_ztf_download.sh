#!/bin/bash
#
# Launch ZTF lightcurve download for SDSS DR16Q × Gaia sources
#
# Downloads ZTF lightcurves for 223,219 sources with Gaia epoch photometry
# Runtime: Several hours to days depending on API performance
#

set -e

# Activate virtual environment
source ~/Work/venvs/.venv/bin/activate

# Create output directory
# mkdir -p ~/data/ztf

# Create local chunks directory
# mkdir -p ztf_chunks

# Launch download
python download_ztf_threaded.py \
  --input data/sdssdr16q_gaia_ztf_targets.fits \
  --output ~/data/ztf/sdssdr16q_ztf_lightcurves.fits \
  --chunks-dir ztf_chunks \
  --manifest ztf_manifest.json \
  --log ztf_download.log \
  --batch-size 500 \
  --threads 8

echo ""
echo "ZTF download complete!"
echo "Output: ~/data/ztf/sdssdr16q_ztf_lightcurves.fits"
