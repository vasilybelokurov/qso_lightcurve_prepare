#!/usr/bin/env python
"""
Clean ZTF lightcurves by applying quality cuts.

Applies quality cuts to g and r bands:
- catflags == 0 (no quality flags)
- |sharp| < 0.25 (PSF sharpness)
- 0.5 <= chi <= 1.5 (PSF fit quality)
- airmass < 1.8 (observation quality)

Based on quality cuts from qso_gp_mock/ztf_10k_lightcurve_quality.ipynb

Usage:
    python clean_ztf_lightcurves.py \
        --input data/gaia_ztf_qso_sample_ztf_lc.fits \
        --output data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits \
        --bands g,r
"""

import os
import sys
import numpy as np
from astropy.table import Table
from tqdm import tqdm
import argparse


def apply_quality_cuts(ztf_table, bands=['g', 'r']):
    """
    Apply quality cuts to ZTF lightcurves for specified bands.

    Quality cuts (per band):
    - catflags == 0
    - |sharp| < 0.25
    - 0.5 <= chi <= 1.5
    - airmass < 1.8

    Parameters
    ----------
    ztf_table : astropy.table.Table
        Original ZTF lightcurve table
    bands : list of str
        Bands to clean ('g', 'r')

    Returns
    -------
    astropy.table.Table
        New table with cleaned lightcurves (only specified bands)
    """
    n_sources = len(ztf_table)

    print(f"\nApplying quality cuts to {n_sources} sources...")
    print("Quality cuts:")
    print("  - catflags == 0")
    print("  - |sharp| < 0.25")
    print("  - 0.5 <= chi <= 1.5")
    print("  - airmass < 1.8")
    print(f"  - Bands: {', '.join(bands)}")

    # Create new table with metadata columns
    cleaned = Table()
    cleaned['source_id'] = ztf_table['source_id']
    cleaned['ra'] = ztf_table['ra']
    cleaned['dec'] = ztf_table['dec']

    # Process each band
    for band in bands:
        print(f"\nProcessing {band}-band...")

        mjd_clean = []
        mag_clean = []
        magerr_clean = []
        catflags_clean = []
        sharp_clean = []
        chi_clean = []
        limitmag_clean = []
        airmass_clean = []

        n_epochs_before = []
        n_epochs_after = []

        for i in tqdm(range(n_sources), desc=f"Cleaning {band}", ncols=80):
            # Apply quality cuts
            good = (ztf_table[f'catflags_{band}'][i] == 0)
            good &= (np.abs(ztf_table[f'sharp_{band}'][i]) < 0.25)
            good &= (ztf_table[f'chi_{band}'][i] >= 0.5)
            good &= (ztf_table[f'chi_{band}'][i] <= 1.5)
            good &= (ztf_table[f'airmass_{band}'][i] < 1.8)

            n_epochs_before.append(len(ztf_table[f'mjd_{band}'][i]))
            n_epochs_after.append(good.sum())

            # Keep cleaned arrays
            mjd_clean.append(np.array(ztf_table[f'mjd_{band}'][i][good], dtype=float))
            mag_clean.append(np.array(ztf_table[f'mag_{band}'][i][good], dtype=np.float32))
            magerr_clean.append(np.array(ztf_table[f'magerr_{band}'][i][good], dtype=np.float32))
            catflags_clean.append(np.array(ztf_table[f'catflags_{band}'][i][good], dtype=np.int32))
            sharp_clean.append(np.array(ztf_table[f'sharp_{band}'][i][good], dtype=np.float32))
            chi_clean.append(np.array(ztf_table[f'chi_{band}'][i][good], dtype=np.float32))
            limitmag_clean.append(np.array(ztf_table[f'limitmag_{band}'][i][good], dtype=np.float32))
            airmass_clean.append(np.array(ztf_table[f'airmass_{band}'][i][good], dtype=np.float32))

        # Add to cleaned table
        cleaned[f'mjd_{band}'] = mjd_clean
        cleaned[f'mag_{band}'] = mag_clean
        cleaned[f'magerr_{band}'] = magerr_clean
        cleaned[f'catflags_{band}'] = catflags_clean
        cleaned[f'sharp_{band}'] = sharp_clean
        cleaned[f'chi_{band}'] = chi_clean
        cleaned[f'limitmag_{band}'] = limitmag_clean
        cleaned[f'airmass_{band}'] = airmass_clean

        # Statistics
        n_epochs_before_arr = np.array(n_epochs_before)
        n_epochs_after_arr = np.array(n_epochs_after)

        mean_before = np.mean(n_epochs_before_arr)
        mean_after = np.mean(n_epochs_after_arr)
        median_before = np.median(n_epochs_before_arr)
        median_after = np.median(n_epochs_after_arr)

        n_zero = np.sum(n_epochs_after_arr == 0)
        retention_rate = mean_after / mean_before * 100 if mean_before > 0 else 0

        print(f"  {band}-band statistics:")
        print(f"    Mean epochs: {mean_before:.1f} → {mean_after:.1f} ({retention_rate:.1f}% retained)")
        print(f"    Median epochs: {median_before:.0f} → {median_after:.0f}")
        print(f"    Sources with 0 epochs after cleaning: {n_zero}/{n_sources} ({100*n_zero/n_sources:.1f}%)")

    return cleaned


def main():
    parser = argparse.ArgumentParser(
        description="Clean ZTF lightcurves with quality cuts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input',
        type=str,
        default='data/gaia_ztf_qso_sample_ztf_lc.fits',
        help='Input ZTF lightcurve FITS file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits',
        help='Output cleaned FITS file'
    )
    parser.add_argument(
        '--bands',
        type=str,
        default='g,r',
        help='Comma-separated bands to clean'
    )

    args = parser.parse_args()

    # Parse bands
    bands = args.bands.split(',')
    bands = [b.strip().lower() for b in bands]

    print("="*70)
    print("ZTF LIGHTCURVE QUALITY CLEANING")
    print("="*70)
    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Bands:  {', '.join(bands)}")

    # Load ZTF catalog
    print(f"\nLoading ZTF catalog from {args.input}...")
    ztf = Table.read(args.input)
    print(f"  Loaded {len(ztf)} sources")
    print(f"  Columns: {len(ztf.colnames)}")

    # Apply quality cuts
    cleaned = apply_quality_cuts(ztf, bands=bands)

    # Save cleaned catalog
    print(f"\n{'='*70}")
    print("SAVING CLEANED CATALOG")
    print("="*70)
    print(f"Output file: {args.output}")
    print(f"Sources: {len(cleaned)}")
    print(f"Columns: {len(cleaned.colnames)} ({', '.join(cleaned.colnames[:5])}...)")

    cleaned.write(args.output, format='fits', overwrite=True)

    file_size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"✓ Saved {len(cleaned)} sources ({file_size_mb:.1f} MB)")

    print("\n" + "="*70)
    print("CLEANING COMPLETE")
    print("="*70)
    print(f"Original file preserved: {args.input}")
    print(f"Cleaned file created: {args.output}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
