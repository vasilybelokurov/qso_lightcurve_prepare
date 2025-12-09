#!/usr/bin/env python3
"""
Download Gaia DR3 epoch photometry for SDSS DR16Q × Gaia cross-matched sources.

Uses local_join to upload source_ids from Step 1 catalog and retrieve
epoch photometry from WSDB.

Requires: ~/data/gaia/sdssdr16q_gaia_source.fits (from Step 1)
Output: ~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits
"""

import sqlutilpy as sqlutil
from astropy.table import Table
import os
import time


def download_epoch_photometry():
    """Download epoch photometry using local_join."""

    print("=" * 70)
    print("SDSS DR16Q × Gaia DR3 Epoch Photometry Download")
    print("=" * 70)
    print()

    # Load Step 1 catalog
    catalog_file = os.path.expanduser('~/data/gaia/sdssdr16q_gaia_source.fits')

    if not os.path.exists(catalog_file):
        print(f"ERROR: Catalog file not found: {catalog_file}")
        print("Please run download_sdssdr16q_gaia_source.py first.")
        return 1

    print(f"Step 1: Loading catalog from {catalog_file}...")
    catalog = Table.read(catalog_file)
    print(f"  ✓ Loaded {len(catalog)} sources")
    print()

    # Extract source_ids
    if 'source_id' not in catalog.colnames:
        print("ERROR: 'source_id' column not found in catalog")
        return 1

    source_ids = catalog['source_id']
    print(f"Step 2: Preparing {len(source_ids)} source_ids for upload...")
    print()

    # Use local_join (like Quaia pattern)
    print("Step 3: Uploading source_ids and downloading epoch photometry...")
    print("  This will take 10-30 minutes depending on server load...")
    print()

    t0 = time.time()

    result = sqlutil.local_join(
        '''
        SELECT
            m.source_id,
            e.g_transit_time, e.g_transit_mag, e.g_transit_flux, e.g_transit_flux_error,
            e.bp_obs_time, e.bp_mag, e.bp_flux, e.bp_flux_error,
            e.rp_obs_time, e.rp_mag, e.rp_flux, e.rp_flux_error,
            e.photometry_flag_noisy_data,
            e.photometry_flag_sm_unavailable,
            e.photometry_flag_bp_unavailable,
            e.photometry_flag_rp_unavailable,
            e.photometry_flag_sm_reject,
            e.variability_flag_g_reject,
            e.variability_flag_bp_reject,
            e.variability_flag_rp_reject
        FROM mytmptable AS m
        JOIN gaia_dr3.epoch_photometry AS e ON e.source_id = m.source_id
        ''',
        'mytmptable',
        (source_ids,),
        ('source_id',),
        asDict=True
    )

    elapsed = time.time() - t0

    print(f"  ✓ Downloaded epoch photometry in {elapsed/60:.1f} minutes")
    print()

    # Convert to table
    epoch_table = Table(result)

    print(f"Step 4: Saving epoch photometry...")
    output_file = os.path.expanduser('~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits')

    epoch_table.write(output_file, overwrite=True)

    # Get file size
    file_size_mb = os.path.getsize(output_file) / (1024 * 1024)

    print(f"  ✓ Saved: {output_file}")
    print(f"  File size: {file_size_mb:.1f} MB")
    print(f"  Rows: {len(epoch_table)}")
    print(f"  Columns: {len(epoch_table.colnames)}")
    print()

    # Print summary statistics
    print("=" * 70)
    print("Summary Statistics")
    print("=" * 70)

    unique_sources = len(set(epoch_table['source_id']))
    print(f"Total epoch photometry records: {len(epoch_table)}")
    print(f"Unique sources with epochs: {unique_sources}")
    print(f"Average epochs per source: {len(epoch_table)/unique_sources:.1f}")

    # Count epochs per band
    import numpy as np

    def count_band_epochs(times_col):
        """Count number of epochs per source."""
        return [len(t) if t is not None else 0 for t in epoch_table[times_col]]

    if 'g_transit_time' in epoch_table.colnames:
        g_epochs = count_band_epochs('g_transit_time')
        print(f"G-band: median={np.median(g_epochs):.0f} epochs, total={np.sum(g_epochs)}")

    if 'bp_obs_time' in epoch_table.colnames:
        bp_epochs = count_band_epochs('bp_obs_time')
        print(f"BP-band: median={np.median(bp_epochs):.0f} epochs, total={np.sum(bp_epochs)}")

    if 'rp_obs_time' in epoch_table.colnames:
        rp_epochs = count_band_epochs('rp_obs_time')
        print(f"RP-band: median={np.median(rp_epochs):.0f} epochs, total={np.sum(rp_epochs)}")

    print("=" * 70)
    print()
    print("Done! You now have:")
    print(f"  1. {catalog_file}")
    print(f"  2. {output_file}")

    return 0


def main():
    return download_epoch_photometry()


if __name__ == '__main__':
    import sys
    sys.exit(main())
