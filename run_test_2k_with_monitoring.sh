#!/bin/bash
#
# Test ZTF download on first 2K sources with performance monitoring
#
# Monitors:
# - Download speed (sources/sec)
# - CPU usage
# - Network activity
# - Memory usage
#

set -e

# Activate virtual environment
source ~/Work/venvs/.venv/bin/activate

# Create test catalog if needed
if [ ! -f data/test_2k_ztf_targets.fits ]; then
    echo "Creating test catalog (first 2000 sources)..."
    python create_test_2k_catalog.py
fi

# Create output directories
mkdir -p test_2k_chunks

# Clean up previous test run
rm -f test_2k_manifest.json
rm -rf test_2k_chunks/*

echo ""
echo "========================================================================"
echo "ZTF Download Test: 2000 Sources"
echo "========================================================================"
echo ""
echo "Configuration:"
echo "  Input: data/test_2k_ztf_targets.fits (2000 sources)"
echo "  Batch size: 100"
echo "  Threads: 8"
echo "  Output: test_2k_chunks/"
echo ""
echo "Monitoring CPU usage every 5 seconds..."
echo "Press Ctrl+C to stop monitoring (download will continue)"
echo ""
echo "========================================================================"
echo ""

# Get Python PID for monitoring
DOWNLOAD_LOG="test_2k_download.log"

# Start download in background
python download_ztf_threaded.py \
  --input data/test_2k_ztf_targets.fits \
  --output test_2k_ztf_lightcurves.fits \
  --chunks-dir test_2k_chunks \
  --manifest test_2k_manifest.json \
  --log "$DOWNLOAD_LOG" \
  --batch-size 100 \
  --threads 8 &

DOWNLOAD_PID=$!

echo "Download started (PID: $DOWNLOAD_PID)"
echo ""

# Monitor CPU and progress
monitor_stats() {
    echo "Time | CPU% | Matched | Rate | Status"
    echo "-----+------+---------+------+------------------"

    start_time=$(date +%s)

    while kill -0 $DOWNLOAD_PID 2>/dev/null; do
        current_time=$(date +%s)
        elapsed=$((current_time - start_time))

        # Get CPU usage for Python process
        cpu=$(ps -p $DOWNLOAD_PID -o %cpu= 2>/dev/null | xargs || echo "0.0")

        # Parse log for progress
        if [ -f "$DOWNLOAD_LOG" ]; then
            # Get last "Total: X/2000" line
            matched=$(grep -o "Total: [0-9]*/2000" "$DOWNLOAD_LOG" 2>/dev/null | tail -1 | cut -d' ' -f2 | cut -d'/' -f1 || echo "0")

            # Calculate rate
            if [ "$matched" != "0" ] && [ "$elapsed" -gt 0 ]; then
                rate=$(echo "scale=2; $matched / $elapsed" | bc)
            else
                rate="0.00"
            fi

            # Get last batch status
            status=$(tail -1 "$DOWNLOAD_LOG" 2>/dev/null | cut -c1-40 || echo "Starting...")
        else
            matched="0"
            rate="0.00"
            status="Initializing..."
        fi

        # Format time
        mins=$((elapsed / 60))
        secs=$((elapsed % 60))
        time_str=$(printf "%02d:%02d" $mins $secs)

        printf "%5s | %5s | %7s | %4s | %-40s\n" "$time_str" "$cpu" "$matched" "$rate" "$status"

        sleep 5
    done
}

# Run monitoring
monitor_stats

# Wait for download to complete
wait $DOWNLOAD_PID
EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Download Complete"
echo "========================================================================"
echo ""

# Show final statistics
if [ -f "$DOWNLOAD_LOG" ]; then
    echo "Final statistics from log:"
    tail -5 "$DOWNLOAD_LOG"
    echo ""
fi

# Show chunk files
if [ -d test_2k_chunks ]; then
    echo "Chunk files created:"
    ls -lh test_2k_chunks/*.fits 2>/dev/null | head -10
    chunk_count=$(ls test_2k_chunks/*.fits 2>/dev/null | wc -l)
    echo "  Total chunks: $chunk_count"
    echo ""
fi

# Show manifest
if [ -f test_2k_manifest.json ]; then
    echo "Manifest summary:"
    cat test_2k_manifest.json
    echo ""
fi

echo "Test log: $DOWNLOAD_LOG"
echo ""

exit $EXIT_CODE
