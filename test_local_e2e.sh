#!/bin/bash
# Local end-to-end testing script (synced with GitHub Actions workflow)
# Usage: bash test_local_e2e.sh /path/to/mypipelinerun_circinus_p3_3-MFS-image.fits
#
# NOTE: Requires catalog dependencies:
#   pip install breizorro[catalog]  # For cataloging features
#   pip install breizorro[all]      # For all features (catalog + GUI)

FITS_FILE="${1:-.}/mypipelinerun_circinus_p3_3-MFS-image.fits"

if [ ! -f "$FITS_FILE" ]; then
    echo "Downloading test dataset..."
    curl -L -o mypipelinerun_circinus_p3_3-MFS-image.tar \
        "https://www.dropbox.com/scl/fi/xoaqth4jhjnx098xbegvn/mypipelinerun_circinus_p3_3-MFS-image.tar?rlkey=vdgkz2xpjgpp7v7hjxm7kg6g2&st=isg70eve&dl=1"
    tar -xvf mypipelinerun_circinus_p3_3-MFS-image.tar
fi

echo "========== Basic mask generation =========="
uv run breizorro -r "$FITS_FILE" --outfile test.mask.fits --outregion test.mask.rgn

echo "========== Test with remove islands and numbering =========="
uv run breizorro -r "$FITS_FILE" -t 10 --remove-islands 19,21 --number-islands --outfile test.numbered.fits

echo "========== Test catalog generation with centroid fitting (default) =========="
uv run breizorro -r "$FITS_FILE" --outcatalog test_centroid.cat --source-fitting centroid --outfile test_centroid.fits 2>&1 | grep -i "welcome\|catalog\|error" || true

echo "========== Test catalog generation with Gaussian fitting =========="
uv run breizorro -r "$FITS_FILE" --outcatalog test_gaussian.cat --source-fitting gaussian --outfile test_gaussian.fits 2>&1 | grep -i "welcome\|catalog\|error" || true

echo "========== Test radial cutoff feature (beam attenuation) =========="
uv run breizorro -r "$FITS_FILE" --radial-cutoff 400 --outfile test.radial.fits

echo "========== Test dilation =========="
uv run breizorro -r "$FITS_FILE" --dilate 3 --outfile test.dilate.fits

echo ""
echo "✅ All basic tests completed!"

# Cleanup test artifacts
echo ""
echo "Cleaning up test artifacts..."
rm -f test*.fits test*.cat test*.rgn
echo "✅ Cleanup done"
