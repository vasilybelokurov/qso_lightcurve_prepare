#!/usr/bin/env python
"""
Create matched ZTF+Gaia+Catalog files.

Takes ZTF lightcurve data and creates corresponding Gaia lightcurves
and catalog files for the same sources (matched by source_id).

Adapted from qso_gp_mock/create_matched_ztf_gaia_catalog.py

Usage:
    python create_matched_catalog.py \
        --ztf-lc data/ztf_lightcurves_partial_merged.fits \
        --gaia-lc ~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits \
        --gaia-catalog ~/data/gaia/sdssdr16q_gaia_source.fits \
        --output-prefix data/gaia_ztf_qso_sample
"""

import os
import sys
import numpy as np
from astropy.table import Table
import argparse


def load_ztf_lightcurves(ztf_path):
    """Load ZTF lightcurves."""
    print(f"Loading ZTF lightcurves from {ztf_path}...")
    ztf = Table.read(os.path.expanduser(ztf_path))
    print(f"  Loaded {len(ztf)} ZTF sources")
    return ztf


def load_gaia_lightcurves(gaia_path, source_ids):
    """
    Load Gaia lightcurves for specific source_ids.

    Parameters
    ----------
    gaia_path : str
        Path to Gaia lightcurve FITS file
    source_ids : array-like
        List of Gaia DR3 source_ids to extract

    Returns
    -------
    gaia_lc : astropy.table.Table
        Gaia lightcurves for matched sources
    """
    print(f"Loading Gaia lightcurves from {gaia_path}...")

    # Load full Gaia catalog as Table (handles variable-length arrays correctly)
    gaia_full = Table.read(os.path.expanduser(gaia_path))

    print(f"  Total Gaia sources: {len(gaia_full)}")
    print(f"  Looking for {len(source_ids)} matches...")

    # Create index for fast lookup
    source_id_set = set(source_ids)
    mask = np.array([sid in source_id_set for sid in gaia_full['source_id']])

    print(f"  Found {mask.sum()} matches")

    # Extract matched rows
    gaia_lc = gaia_full[mask]

    # Sort by source_id to match ZTF order
    gaia_lc.sort('source_id')

    return gaia_lc


def load_gaia_catalog(catalog_path, source_ids):
    """
    Load Gaia source catalog data for specific source_ids.

    Parameters
    ----------
    catalog_path : str
        Path to Gaia catalog FITS file
    source_ids : array-like
        List of Gaia DR3 source_ids to extract

    Returns
    -------
    catalog : astropy.table.Table
        Catalog data for matched sources (ALL columns)
    """
    print(f"Loading Gaia source catalog from {catalog_path}...")

    # Load catalog
    catalog_full = Table.read(os.path.expanduser(catalog_path))
    print(f"  Total catalog sources: {len(catalog_full)}")
    print(f"  Total columns: {len(catalog_full.colnames)}")

    # Find matches
    source_id_set = set(source_ids)
    mask = np.array([sid in source_id_set for sid in catalog_full['source_id']])

    print(f"  Found {mask.sum()} matches")

    catalog = catalog_full[mask]

    # Sort by source_id to match ZTF order
    catalog.sort('source_id')

    return catalog


def create_fully_matched_subset(ztf, gaia_lc, catalog):
    """
    Create subset where all three datasets have data.

    Returns
    -------
    ztf_matched, gaia_matched, catalog_matched : Tables
        Tables with only sources that exist in all three datasets
    """
    print("\nCreating fully-matched subset...")

    # Find common source_ids
    ztf_ids = set(ztf['source_id'])
    gaia_ids = set(gaia_lc['source_id'])
    catalog_ids = set(catalog['source_id'])

    common_ids = ztf_ids & gaia_ids & catalog_ids
    print(f"  ZTF sources: {len(ztf_ids)}")
    print(f"  Gaia LC sources: {len(gaia_ids)}")
    print(f"  Catalog sources: {len(catalog_ids)}")
    print(f"  Common sources (ZTF ∩ Gaia LC ∩ Catalog): {len(common_ids)}")

    if len(common_ids) == 0:
        print("  ERROR: No sources with all three datasets!")
        return None, None, None

    # Filter each table
    ztf_mask = np.array([sid in common_ids for sid in ztf['source_id']])
    gaia_mask = np.array([sid in common_ids for sid in gaia_lc['source_id']])
    catalog_mask = np.array([sid in common_ids for sid in catalog['source_id']])

    ztf_matched = ztf[ztf_mask]
    gaia_matched = gaia_lc[gaia_mask]
    catalog_matched = catalog[catalog_mask]

    # Sort all by source_id
    ztf_matched.sort('source_id')
    gaia_matched.sort('source_id')
    catalog_matched.sort('source_id')

    # Verify alignment
    assert len(ztf_matched) == len(gaia_matched) == len(catalog_matched), \
        "Row counts don't match after filtering!"
    assert np.all(ztf_matched['source_id'] == gaia_matched['source_id']), \
        "ZTF and Gaia source_ids not aligned!"
    assert np.all(ztf_matched['source_id'] == catalog_matched['source_id']), \
        "ZTF and Catalog source_ids not aligned!"

    print(f"  ✓ Fully-matched subset: {len(ztf_matched)} sources aligned")
    print(f"  ✓ All three tables have identical row counts")
    print(f"  ✓ source_id alignment verified")

    return ztf_matched, gaia_matched, catalog_matched


def print_summary(ztf, gaia_lc, catalog):
    """Print summary statistics."""
    print("\n" + "="*70)
    print("MATCHED CATALOG SUMMARY")
    print("="*70)

    print(f"\nNumber of matched sources: {len(ztf)}")

    # ZTF statistics
    print(f"\nZTF lightcurves:")
    n_ztf_g = [len(ztf['mjd_g'][i]) for i in range(len(ztf))]
    n_ztf_r = [len(ztf['mjd_r'][i]) for i in range(len(ztf))]
    n_ztf_i = [len(ztf['mjd_i'][i]) for i in range(len(ztf))]

    n_with_g = sum(1 for n in n_ztf_g if n > 0)
    n_with_r = sum(1 for n in n_ztf_r if n > 0)
    n_with_i = sum(1 for n in n_ztf_i if n > 0)

    print(f"  g-band: {n_with_g} sources, median {np.median(n_ztf_g):.0f} epochs (max {np.max(n_ztf_g)})")
    print(f"  r-band: {n_with_r} sources, median {np.median(n_ztf_r):.0f} epochs (max {np.max(n_ztf_r)})")
    print(f"  i-band: {n_with_i} sources, median {np.median(n_ztf_i):.0f} epochs (max {np.max(n_ztf_i)})")
    print(f"  Total ZTF epochs: {np.sum(n_ztf_g) + np.sum(n_ztf_r) + np.sum(n_ztf_i):,}")

    # Gaia statistics
    print(f"\nGaia lightcurves:")
    n_gaia_g = [len(gaia_lc['g_transit_time'][i]) for i in range(len(gaia_lc))]
    n_gaia_bp = [len(gaia_lc['bp_obs_time'][i]) for i in range(len(gaia_lc))]
    n_gaia_rp = [len(gaia_lc['rp_obs_time'][i]) for i in range(len(gaia_lc))]

    print(f"  g-band: median {np.median(n_gaia_g):.0f} epochs (max {np.max(n_gaia_g)})")
    print(f"  bp-band: median {np.median(n_gaia_bp):.0f} epochs (max {np.max(n_gaia_bp)})")
    print(f"  rp-band: median {np.median(n_gaia_rp):.0f} epochs (max {np.max(n_gaia_rp)})")
    print(f"  Total Gaia epochs: {np.sum(n_gaia_g) + np.sum(n_gaia_bp) + np.sum(n_gaia_rp):,}")

    # Catalog statistics
    print(f"\nCatalog:")
    print(f"  Columns: {len(catalog.colnames)}")
    if 'z' in catalog.colnames:
        print(f"  Redshift range: z={catalog['z'].min():.3f} - {catalog['z'].max():.3f}")
        print(f"  Median redshift: z={np.median(catalog['z']):.3f}")

    print("="*70)


def save_files(ztf, gaia_lc, catalog, output_prefix):
    """Save matched files."""

    # ZTF lightcurves
    ztf_output = f"{output_prefix}_ztf_lc.fits"
    print(f"\nSaving ZTF lightcurves to {ztf_output}...")
    ztf.write(ztf_output, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(ztf_output) / (1024 * 1024)
    print(f"  ✓ Saved {len(ztf)} sources ({file_size_mb:.1f} MB)")

    # Gaia lightcurves
    gaia_output = f"{output_prefix}_gaia_lc.fits"
    print(f"\nSaving Gaia lightcurves to {gaia_output}...")
    gaia_lc.write(gaia_output, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(gaia_output) / (1024 * 1024)
    print(f"  ✓ Saved {len(gaia_lc)} sources ({file_size_mb:.1f} MB)")

    # Catalog
    catalog_output = f"{output_prefix}_catalog.fits"
    print(f"\nSaving catalog to {catalog_output}...")
    catalog.write(catalog_output, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(catalog_output) / (1024 * 1024)
    print(f"  ✓ Saved {len(catalog)} sources, {len(catalog.colnames)} columns ({file_size_mb:.1f} MB)")

    print(f"\n" + "="*70)
    print("FILES CREATED:")
    print(f"  {ztf_output}")
    print(f"  {gaia_output}")
    print(f"  {catalog_output}")
    print("="*70)

    return ztf_output, gaia_output, catalog_output


def main():
    parser = argparse.ArgumentParser(
        description="Create matched ZTF+Gaia+Catalog files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--ztf-lc',
        type=str,
        default='data/ztf_lightcurves_partial_merged.fits',
        help='Input ZTF lightcurves FITS file'
    )
    parser.add_argument(
        '--gaia-lc',
        type=str,
        default='~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits',
        help='Gaia epoch photometry FITS file'
    )
    parser.add_argument(
        '--gaia-catalog',
        type=str,
        default='~/data/gaia/sdssdr16q_gaia_source.fits',
        help='Gaia source catalog FITS file'
    )
    parser.add_argument(
        '--output-prefix',
        type=str,
        default='data/gaia_ztf_qso_sample',
        help='Output filename prefix'
    )

    args = parser.parse_args()

    print("="*70)
    print("CREATE MATCHED ZTF+GAIA+CATALOG FILES")
    print("="*70)

    # Load ZTF data (this defines which sources we have)
    ztf = load_ztf_lightcurves(args.ztf_lc)
    source_ids = ztf['source_id']

    # Load matching Gaia lightcurves
    gaia_lc = load_gaia_lightcurves(args.gaia_lc, source_ids)

    # Load matching catalog
    catalog = load_gaia_catalog(args.gaia_catalog, source_ids)

    # Create fully-matched subset (all three datasets)
    ztf_matched, gaia_matched, catalog_matched = create_fully_matched_subset(
        ztf, gaia_lc, catalog
    )

    if ztf_matched is None:
        print("\nERROR: No common sources found!")
        return 1

    # Print summary
    print_summary(ztf_matched, gaia_matched, catalog_matched)

    # Save files
    save_files(ztf_matched, gaia_matched, catalog_matched, args.output_prefix)

    print("\n✓ All done!")
    print(f"\nMatching summary:")
    print(f"  Input ZTF sources: {len(ztf)}")
    print(f"  Found in Gaia LC: {len(gaia_lc)}")
    print(f"  Found in Catalog: {len(catalog)}")
    print(f"  Final matched (ALL THREE): {len(ztf_matched)}")
    print(f"  Match rate: {100*len(ztf_matched)/len(ztf):.1f}%")

    return 0


if __name__ == '__main__':
    sys.exit(main())
