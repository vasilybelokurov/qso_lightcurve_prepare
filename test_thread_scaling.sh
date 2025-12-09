#!/bin/bash
#
# Test different thread counts to find optimal configuration
#

set -e

source ~/Work/venvs/.venv/bin/activate

# Create test catalog if needed
if [ ! -f data/test_2k_ztf_targets.fits ]; then
    python create_test_2k_catalog.py
fi

# Test with 100 sources at different thread counts
test_threads() {
    local n_threads=$1
    local test_size=100

    echo ""
    echo "Testing with $n_threads threads ($test_size sources)..."

    # Clean up
    rm -rf test_threads_${n_threads}
    mkdir -p test_threads_${n_threads}

    # Extract first 100 sources
    python -c "
from astropy.table import Table
cat = Table.read('data/test_2k_ztf_targets.fits')
cat[:${test_size}].write('test_threads_${n_threads}/input.fits', overwrite=True)
"

    # Time the download
    start=$(date +%s)

    python download_ztf_threaded.py \
      --input test_threads_${n_threads}/input.fits \
      --output test_threads_${n_threads}/output.fits \
      --chunks-dir test_threads_${n_threads}/chunks \
      --manifest test_threads_${n_threads}/manifest.json \
      --log test_threads_${n_threads}/download.log \
      --batch-size ${test_size} \
      --threads ${n_threads} > /dev/null 2>&1

    end=$(date +%s)
    elapsed=$((end - start))

    # Get match count
    matched=$(grep "Total:" test_threads_${n_threads}/download.log | tail -1 | grep -o "[0-9]*/100" | cut -d'/' -f1)

    if [ -z "$matched" ]; then
        matched=0
    fi

    rate=$(echo "scale=2; $matched / $elapsed" | bc)

    echo "  Threads: $n_threads | Time: ${elapsed}s | Matched: $matched/$test_size | Rate: $rate src/s"
}

echo "========================================"
echo "Thread Scaling Test (100 sources each)"
echo "========================================"

# Test different thread counts
for threads in 4 8 16 24 32; do
    test_threads $threads
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo "Compare the rates above to find optimal thread count"
