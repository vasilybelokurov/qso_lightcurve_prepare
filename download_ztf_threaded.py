#!/usr/bin/env python
"""
Threaded ZTF lightcurve downloader (faster for I/O-bound API calls).

Features:
- Multi-threaded downloading (better for I/O-bound operations)
- Saves each batch as separate FITS file
- Manifest file tracks completed batches
- Resume capability
- Thread-safe logging and file writing
"""

import os
import sys
import time
import json
import numpy as np
import requests
from io import StringIO
from astropy.table import Table, vstack
from astropy.io import fits, ascii
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import warnings
warnings.filterwarnings('ignore')


# IRSA ZTF API endpoint
ZTF_API_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"

# Configuration
SEARCH_RADIUS_ARCSEC = 2.0
RETRY_ATTEMPTS = 3
RETRY_DELAY = 2.0
QUERY_TIMEOUT = 30
BATCH_SIZE = 100
N_THREADS = 8  # More threads for I/O-bound operations

# Thread-safe lock for logging
log_lock = threading.Lock()


def setup_logging(log_file):
    """Setup logging to both file and console."""
    logger = logging.getLogger('ztf_threaded')
    logger.setLevel(logging.INFO)
    logger.handlers = []

    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def load_catalog(fits_path):
    """Load Gaia quasar catalog."""
    cat = Table.read(os.path.expanduser(fits_path))
    required_cols = ['source_id', 'ra', 'dec']
    for col in required_cols:
        if col not in cat.colnames:
            raise ValueError(f"Missing required column: {col}")
    return cat


def load_manifest(manifest_path):
    """Load manifest of completed batches."""
    if not os.path.exists(manifest_path):
        return {'completed_batches': [], 'n_matched': 0, 'n_failed': 0}

    with open(manifest_path, 'r') as f:
        return json.load(f)


def save_manifest(manifest_path, manifest):
    """Save manifest."""
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)


def query_ztf_lightcurve(ra, dec, radius_arcsec=2.0, bandname='g,r,i', timeout=30):
    """Query ZTF API for single source (copied from working sequential version)."""
    radius_deg = radius_arcsec / 3600.0

    params = {
        'POS': f'CIRCLE {ra:.6f} {dec:.6f} {radius_deg:.6f}',
        'BANDNAME': bandname,
        'FORMAT': 'CSV',
    }

    try:
        response = requests.get(ZTF_API_URL, params=params, timeout=timeout)
        response.raise_for_status()

        csv_text = response.text

        if len(csv_text.strip()) == 0 or 'No rows returned' in csv_text:
            return None

        table = ascii.read(csv_text, format='csv')

        if len(table) == 0:
            return None

        return table

    except requests.exceptions.RequestException:
        return None


def process_source(source):
    """
    Process a single source.

    Parameters
    ----------
    source : astropy.table.Row
        Source row with source_id, ra, dec

    Returns
    -------
    dict or None
        Result dictionary or None if failed
    """
    source_id = int(source['source_id'])
    ra = float(source['ra'])
    dec = float(source['dec'])

    # Try with retries
    lc_data = None
    for attempt in range(RETRY_ATTEMPTS):
        lc_data = query_ztf_lightcurve(ra, dec, timeout=QUERY_TIMEOUT)
        if lc_data is not None:
            break
        time.sleep(RETRY_DELAY)

    if lc_data is None:
        return None

    # Aggregate by filter
    result_row = {
        'source_id': source_id,
        'ra': ra,
        'dec': dec,
    }

    for band in ['zg', 'zr', 'zi']:
        band_short = band[1]  # 'g', 'r', 'i'
        mask = lc_data['filtercode'] == band

        if mask.sum() == 0:
            result_row[f'mjd_{band_short}'] = np.array([], dtype=float)
            result_row[f'mag_{band_short}'] = np.array([], dtype=np.float32)
            result_row[f'magerr_{band_short}'] = np.array([], dtype=np.float32)
            result_row[f'catflags_{band_short}'] = np.array([], dtype=np.int32)
            result_row[f'sharp_{band_short}'] = np.array([], dtype=np.float32)
            result_row[f'chi_{band_short}'] = np.array([], dtype=np.float32)
            result_row[f'limitmag_{band_short}'] = np.array([], dtype=np.float32)
            result_row[f'airmass_{band_short}'] = np.array([], dtype=np.float32)
        else:
            result_row[f'mjd_{band_short}'] = np.array(lc_data['mjd'][mask], dtype=float)
            result_row[f'mag_{band_short}'] = np.array(lc_data['mag'][mask], dtype=np.float32)
            result_row[f'magerr_{band_short}'] = np.array(lc_data['magerr'][mask], dtype=np.float32)
            result_row[f'catflags_{band_short}'] = np.array(lc_data['catflags'][mask], dtype=np.int32)
            result_row[f'sharp_{band_short}'] = np.array(lc_data['sharp'][mask], dtype=np.float32)
            result_row[f'chi_{band_short}'] = np.array(lc_data['chi'][mask], dtype=np.float32)
            result_row[f'limitmag_{band_short}'] = np.array(lc_data['limitmag'][mask], dtype=np.float32)
            result_row[f'airmass_{band_short}'] = np.array(lc_data['airmass'][mask], dtype=np.float32)

    return result_row


def download_batch_threaded(batch_idx, sources, chunks_dir, logger, n_threads=N_THREADS):
    """
    Download a batch of sources using threads.

    Parameters
    ----------
    batch_idx : int
        Batch index
    sources : astropy.table.Table
        Source table
    chunks_dir : str
        Directory for chunks
    logger : logging.Logger
        Logger instance
    n_threads : int
        Number of threads

    Returns
    -------
    dict
        Results dictionary
    """
    n_matched = 0
    n_failed = 0
    batch_results = []

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(process_source, sources[i]): i for i in range(len(sources))}

        for future in as_completed(futures):
            idx = futures[future]
            try:
                result = future.result()

                if result is not None:
                    batch_results.append(result)
                    n_matched += 1
                else:
                    n_failed += 1

                # Log progress every 10 sources
                total_processed = n_matched + n_failed
                if total_processed % 10 == 0:
                    with log_lock:
                        logger.info(f"Batch {batch_idx}: {n_matched}/{total_processed} matched")

            except Exception as e:
                n_failed += 1
                with log_lock:
                    logger.error(f"Batch {batch_idx}: Source {idx} failed: {e}")

    # Save batch chunk
    if len(batch_results) > 0:
        chunk_file = os.path.join(chunks_dir, f'chunk_{batch_idx:04d}.fits')
        batch_table = Table(batch_results)
        batch_table.write(chunk_file, format='fits', overwrite=True)

        with log_lock:
            logger.info(f"Batch {batch_idx}: Saved {n_matched} sources to {chunk_file}")

    return {
        'batch_idx': batch_idx,
        'n_matched': n_matched,
        'n_failed': n_failed
    }


def download_threaded(catalog, chunks_dir, manifest_path, log_file,
                     batch_size=BATCH_SIZE, n_threads=N_THREADS):
    """
    Download ZTF lightcurves using threads.
    """
    logger = setup_logging(log_file)

    # Load manifest
    manifest = load_manifest(manifest_path)
    completed_batches = set(manifest['completed_batches'])

    # Calculate batches
    n_sources = len(catalog)
    n_batches = int(np.ceil(n_sources / batch_size))

    logger.info(f"Starting threaded download: {n_sources} sources, {n_batches} batches, {n_threads} threads")
    logger.info(f"Already completed: {len(completed_batches)} batches")

    # Process batches sequentially (but sources within batch in parallel)
    for batch_idx in range(n_batches):
        if batch_idx in completed_batches:
            logger.info(f"Batch {batch_idx}: Skipping (already completed)")
            continue

        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_sources)
        sources = catalog[start_idx:end_idx]

        logger.info(f"Batch {batch_idx}: Starting ({len(sources)} sources)")
        start_time = time.time()

        result = download_batch_threaded(batch_idx, sources, chunks_dir, logger, n_threads)

        elapsed = time.time() - start_time
        rate = len(sources) / elapsed if elapsed > 0 else 0

        manifest['completed_batches'].append(batch_idx)
        manifest['n_matched'] += result['n_matched']
        manifest['n_failed'] += result['n_failed']
        save_manifest(manifest_path, manifest)

        logger.info(f"Batch {batch_idx}: Complete in {elapsed:.1f}s ({rate:.2f} src/s) | "
                   f"Total: {manifest['n_matched']}/{n_sources}")

    logger.info(f"Download complete: {manifest['n_matched']} matched, {manifest['n_failed']} failed")

    return manifest


def main():
    parser = argparse.ArgumentParser(
        description="Threaded ZTF lightcurve downloader",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', type=str, required=True, help='Input catalog FITS file')
    parser.add_argument('--output', type=str, required=True, help='Output merged FITS file')
    parser.add_argument('--chunks-dir', type=str, required=True, help='Directory for chunk files')
    parser.add_argument('--manifest', type=str, required=True, help='Manifest JSON file')
    parser.add_argument('--log', type=str, required=True, help='Log file')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE, help='Sources per batch')
    parser.add_argument('--threads', type=int, default=N_THREADS, help='Number of threads')

    args = parser.parse_args()

    # Create chunks directory
    os.makedirs(args.chunks_dir, exist_ok=True)

    # Load catalog
    print(f"Loading catalog from {args.input}...")
    catalog = load_catalog(args.input)
    print(f"Loaded {len(catalog)} sources")

    # Download with threads
    manifest = download_threaded(
        catalog,
        args.chunks_dir,
        args.manifest,
        args.log,
        batch_size=args.batch_size,
        n_threads=args.threads
    )

    # Merge chunks if all complete
    expected_batches = int(np.ceil(len(catalog) / args.batch_size))
    if len(manifest['completed_batches']) == expected_batches:
        chunk_files = sorted([
            os.path.join(args.chunks_dir, f)
            for f in os.listdir(args.chunks_dir)
            if f.startswith('chunk_') and f.endswith('.fits')
        ])

        if len(chunk_files) > 0:
            print(f"\nMerging {len(chunk_files)} chunks...")
            tables = [Table.read(f) for f in chunk_files]
            merged = vstack(tables)

            merged.write(args.output, format='fits', overwrite=True)
            print(f"✓ Saved merged catalog: {args.output} ({len(merged)} sources)")
        else:
            print(f"\nNo chunk files found to merge (all sources failed)")
    else:
        print(f"\nDownload incomplete: {len(manifest['completed_batches'])} / {expected_batches} batches done")

    return 0


if __name__ == '__main__':
    sys.exit(main())
