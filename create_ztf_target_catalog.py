#!/usr/bin/env python3
"""
Create compact ZTF target catalog from SDSS DR16Q × Gaia epoch photometry.

Extracts unique source_ids from epoch photometry file and joins with
source catalog to get RA/Dec coordinates.

Input:
  - ~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits
  - ~/data/gaia/sdssdr16q_gaia_source.fits

Output:
  - data/sdssdr16q_gaia_ztf_targets.fits (source_id, ra, dec)
"""

from astropy.table import Table
import numpy as np
import os


def create_ztf_targets():
    """Create compact ZTF target catalog."""

    print("=" * 70)
    print("Creating ZTF Target Catalog")
    print("=" * 70)
    print()

    # Load epoch photometry file
    epoch_file = os.path.expanduser('~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits')
    print(f"Step 1: Loading epoch photometry from {epoch_file}...")
    epoch_phot = Table.read(epoch_file)
    print(f"  ✓ Loaded {len(epoch_phot)} rows")
    print()

    # Get unique source_ids
    print("Step 2: Extracting unique source_ids...")
    unique_source_ids = np.unique(epoch_phot['source_id'])
    print(f"  ✓ Found {len(unique_source_ids)} unique sources")
    print()

    # Load source catalog
    source_file = os.path.expanduser('~/data/gaia/sdssdr16q_gaia_source.fits')
    print(f"Step 3: Loading source catalog from {source_file}...")
    source_cat = Table.read(source_file)
    print(f"  ✓ Loaded {len(source_cat)} sources")
    print()

    # Match source_ids
    print("Step 4: Matching source_ids to get RA/Dec...")

    # Create index mapping for fast lookup
    source_id_to_idx = {sid: i for i, sid in enumerate(source_cat['source_id'])}

    # Build target catalog
    matched_indices = []
    for sid in unique_source_ids:
        if sid in source_id_to_idx:
            matched_indices.append(source_id_to_idx[sid])

    print(f"  ✓ Matched {len(matched_indices)}/{len(unique_source_ids)} sources")
    print()

    # Extract matched rows
    matched_sources = source_cat[matched_indices]

    # Create compact catalog with just source_id, ra, dec
    # Use gaia_ra and gaia_dec (from Gaia, more precise than SDSS)
    target_catalog = Table()
    target_catalog['source_id'] = matched_sources['source_id']
    target_catalog['ra'] = matched_sources['gaia_ra']
    target_catalog['dec'] = matched_sources['gaia_dec']

    # Save
    output_file = 'data/sdssdr16q_gaia_ztf_targets.fits'
    print(f"Step 5: Saving target catalog to {output_file}...")
    os.makedirs('data', exist_ok=True)
    target_catalog.write(output_file, overwrite=True)

    # Get file size
    file_size_kb = os.path.getsize(output_file) / 1024

    print(f"  ✓ Saved: {output_file}")
    print(f"  File size: {file_size_kb:.1f} KB")
    print(f"  Rows: {len(target_catalog)}")
    print(f"  Columns: {len(target_catalog.colnames)}")
    print()

    # Summary statistics
    print("=" * 70)
    print("Summary Statistics")
    print("=" * 70)
    print(f"ZTF target sources: {len(target_catalog)}")
    print(f"RA range: {np.min(target_catalog['ra']):.6f} - {np.max(target_catalog['ra']):.6f}")
    print(f"Dec range: {np.min(target_catalog['dec']):.6f} - {np.max(target_catalog['dec']):.6f}")
    print("=" * 70)
    print()
    print("Next step: Run ZTF lightcurve download")

    return 0


def main():
    return create_ztf_targets()


if __name__ == '__main__':
    import sys
    sys.exit(main())
