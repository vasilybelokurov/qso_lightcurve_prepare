#!/bin/bash
#
# Multi-process ZTF download: Split catalog and run multiple downloaders in parallel
#
# Strategy: Divide catalog into N chunks, run N independent Python processes
# Each process uses 8-16 threads internally
#

set -e

# Configuration
N_PROCESSES=4        # Number of parallel Python processes
THREADS_PER_PROC=8   # Threads per process
INPUT_CATALOG="$1"
OUTPUT_DIR="$2"

if [ -z "$INPUT_CATALOG" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_catalog> <output_dir>"
    echo ""
    echo "Example:"
    echo "  $0 data/sdssdr16q_gaia_ztf_targets.fits ztf_parallel_output"
    echo ""
    echo "Configuration:"
    echo "  Processes: $N_PROCESSES"
    echo "  Threads per process: $THREADS_PER_PROC"
    echo "  Total parallelism: $((N_PROCESSES * THREADS_PER_PROC))"
    exit 1
fi

source ~/Work/venvs/.venv/bin/activate

# Create output structure
mkdir -p "$OUTPUT_DIR"
mkdir -p "${OUTPUT_DIR}/chunks"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/temp_catalogs"

echo "========================================"
echo "Multi-Process ZTF Download"
echo "========================================"
echo "Input: $INPUT_CATALOG"
echo "Output: $OUTPUT_DIR"
echo "Processes: $N_PROCESSES"
echo "Threads/proc: $THREADS_PER_PROC"
echo "Total parallelism: $((N_PROCESSES * THREADS_PER_PROC))"
echo "========================================"
echo ""

# Split catalog into N pieces
echo "Splitting catalog into $N_PROCESSES chunks..."

python -c "
from astropy.table import Table
import numpy as np

cat = Table.read('$INPUT_CATALOG')
n_sources = len(cat)
n_proc = $N_PROCESSES

chunk_size = int(np.ceil(n_sources / n_proc))

for i in range(n_proc):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, n_sources)

    chunk = cat[start:end]
    output = '${OUTPUT_DIR}/temp_catalogs/chunk_{:02d}.fits'.format(i)
    chunk.write(output, overwrite=True)

    print(f'  Chunk {i}: {len(chunk)} sources ({start}-{end}) -> {output}')
"

echo ""
echo "Starting $N_PROCESSES parallel download processes..."
echo ""

# Launch parallel downloaders
PIDS=()
for i in $(seq 0 $((N_PROCESSES - 1))); do
    chunk_file="${OUTPUT_DIR}/temp_catalogs/chunk_$(printf '%02d' $i).fits"

    python download_ztf_threaded.py \
      --input "$chunk_file" \
      --output "${OUTPUT_DIR}/output_$(printf '%02d' $i).fits" \
      --chunks-dir "${OUTPUT_DIR}/chunks" \
      --manifest "${OUTPUT_DIR}/manifest_$(printf '%02d' $i).json" \
      --log "${OUTPUT_DIR}/logs/process_$(printf '%02d' $i).log" \
      --batch-size 500 \
      --threads $THREADS_PER_PROC &

    PIDS+=($!)
    echo "  Process $i: PID ${PIDS[$i]} (log: ${OUTPUT_DIR}/logs/process_$(printf '%02d' $i).log)"
done

echo ""
echo "All processes launched. Monitoring progress..."
echo ""

# Monitor progress
monitor_progress() {
    while true; do
        all_done=true

        echo "$(date '+%H:%M:%S') Progress:"

        for i in $(seq 0 $((N_PROCESSES - 1))); do
            log_file="${OUTPUT_DIR}/logs/process_$(printf '%02d' $i).log"

            if kill -0 ${PIDS[$i]} 2>/dev/null; then
                all_done=false

                # Get last progress line
                if [ -f "$log_file" ]; then
                    progress=$(grep "Total:" "$log_file" 2>/dev/null | tail -1 | grep -o "Total: [0-9]*/[0-9]*" || echo "Starting...")
                else
                    progress="Initializing..."
                fi

                echo "  Process $i: Running - $progress"
            else
                echo "  Process $i: Complete"
            fi
        done

        echo ""

        if [ "$all_done" = true ]; then
            break
        fi

        sleep 10
    done
}

# Run monitoring
monitor_progress

# Wait for all processes
echo "Waiting for all processes to complete..."
for pid in "${PIDS[@]}"; do
    wait $pid
done

echo ""
echo "========================================"
echo "All processes complete!"
echo "========================================"
echo ""

# Show final statistics
echo "Final statistics:"
for i in $(seq 0 $((N_PROCESSES - 1))); do
    log_file="${OUTPUT_DIR}/logs/process_$(printf '%02d' $i).log"
    if [ -f "$log_file" ]; then
        stats=$(grep "Download complete:" "$log_file" | tail -1)
        echo "  Process $i: $stats"
    fi
done

echo ""
echo "Chunk files in: ${OUTPUT_DIR}/chunks/"
echo "Logs in: ${OUTPUT_DIR}/logs/"
echo ""
echo "To merge all chunks into single file, run:"
echo "  python -c \"from astropy.table import Table, vstack; import glob; tables = [Table.read(f) for f in sorted(glob.glob('${OUTPUT_DIR}/chunks/*.fits'))]; merged = vstack(tables); merged.write('${OUTPUT_DIR}/merged_lightcurves.fits', overwrite=True); print(f'Merged {len(merged)} sources')\""
