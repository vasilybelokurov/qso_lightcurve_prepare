#!/usr/bin/env python
"""
Extract ZTF g and r band lightcurves for sources that have Gaia lightcurves.

This creates matching ZTF files for the Gaia lightcurves extracted in
data/gaia_g_for_ztf_primary.fits

Usage:
    python extract_ztf_for_gaia_sources.py
"""

import h5py
import numpy as np
from astropy.table import Table
import argparse
import os
from tqdm import tqdm


def extract_ztf_band(ztf_hdf5_path: str, gaia_sdss_names: set, band: str, output_path: str):
    """
    Extract ZTF lightcurves for sources matching Gaia list.

    Parameters
    ----------
    ztf_hdf5_path : str
        Path to ZTF primary HDF5 file
    gaia_sdss_names : set
        Set of SDSS names from Gaia file
    band : str
        'g' or 'r'
    output_path : str
        Output FITS file path
    """
    print(f"\nExtracting ZTF {band}-band...")
    print(f"  Reading: {ztf_hdf5_path}")

    # Read ZTF HDF5
    with h5py.File(ztf_hdf5_path, 'r') as f:
        n_total = len(f['sdss_name'])
        print(f"  Total sources in ZTF {band}-band: {n_total:,}")

        # Find matching sources
        matched_indices = []
        matched_sdss_names = []

        for i in tqdm(range(n_total), desc=f"  Matching {band}-band"):
            name = f['sdss_name'][i]
            if isinstance(name, bytes):
                name = name.decode().strip()
            else:
                name = str(name).strip()

            if name in gaia_sdss_names:
                matched_indices.append(i)
                matched_sdss_names.append(name)

        n_matched = len(matched_indices)
        print(f"  Matched: {n_matched:,} sources")

        if n_matched == 0:
            print(f"  WARNING: No matches found for {band}-band!")
            return

        # Extract data for matched sources
        print(f"  Extracting data...")

        output_data = {
            'sdss_name': matched_sdss_names,
            'objectid': [],
            'ra': [],
            'dec': [],
            'mjd': [],
            'mag': [],
            'magerr': [],
            'clrcoeff': [],
            'n_epochs_clean': [],
            'n_epochs_raw': [],
            'n_objectids_total': [],
            'median_mag': [],
            'time_baseline': [],
            'separation_arcsec': [],
        }

        for idx in tqdm(matched_indices, desc=f"  Reading {band}-band"):
            output_data['objectid'].append(f['objectid'][idx])
            output_data['ra'].append(f['ra'][idx])
            output_data['dec'].append(f['dec'][idx])
            output_data['mjd'].append(f['mjd'][idx])
            output_data['mag'].append(f['mag'][idx])
            output_data['magerr'].append(f['magerr'][idx])
            output_data['clrcoeff'].append(f['clrcoeff'][idx])
            output_data['n_epochs_clean'].append(f['n_epochs_clean'][idx])
            output_data['n_epochs_raw'].append(f['n_epochs_raw'][idx])
            output_data['n_objectids_total'].append(f['n_objectids_total'][idx])
            output_data['median_mag'].append(f['median_mag'][idx])
            output_data['time_baseline'].append(f['time_baseline'][idx])
            output_data['separation_arcsec'].append(f['separation_arcsec'][idx])

    # Create output table
    print(f"  Creating FITS table...")
    output = Table()
    output['sdss_name'] = output_data['sdss_name']
    output['objectid'] = np.array(output_data['objectid'], dtype=np.int64)
    output['ra'] = np.array(output_data['ra'], dtype=np.float64)
    output['dec'] = np.array(output_data['dec'], dtype=np.float64)
    output['mjd'] = output_data['mjd']
    output['mag'] = output_data['mag']
    output['magerr'] = output_data['magerr']
    output['clrcoeff'] = output_data['clrcoeff']
    output['n_epochs_clean'] = np.array(output_data['n_epochs_clean'], dtype=np.int32)
    output['n_epochs_raw'] = np.array(output_data['n_epochs_raw'], dtype=np.int32)
    output['n_objectids_total'] = np.array(output_data['n_objectids_total'], dtype=np.int32)
    output['median_mag'] = np.array(output_data['median_mag'], dtype=np.float32)
    output['time_baseline'] = np.array(output_data['time_baseline'], dtype=np.float32)
    output['separation_arcsec'] = np.array(output_data['separation_arcsec'], dtype=np.float32)

    # Write to FITS
    print(f"  Writing: {output_path}")
    output.write(output_path, format='fits', overwrite=True)

    file_size_mb = os.path.getsize(output_path) / (1024**2)

    # Statistics
    total_epochs = np.sum(output_data['n_epochs_clean'])
    print(f"\n  {band}-band summary:")
    print(f"    Sources: {n_matched:,}")
    print(f"    Total epochs: {total_epochs:,}")
    print(f"    Median epochs/source: {np.median(output_data['n_epochs_clean']):.0f}")
    print(f"    File size: {file_size_mb:.1f} MB")


def main():
    parser = argparse.ArgumentParser(
        description='Extract ZTF lightcurves for sources with Gaia data'
    )
    parser.add_argument(
        '--gaia-file',
        type=str,
        default='data/gaia_g_for_ztf_primary.fits',
        help='Gaia file with sdss_name list'
    )
    parser.add_argument(
        '--ztf-g',
        type=str,
        default='~/data/ztf/qso/ztf_g_band_lc_primary.h5',
        help='ZTF g-band primary HDF5 file'
    )
    parser.add_argument(
        '--ztf-r',
        type=str,
        default='~/data/ztf/qso/ztf_r_band_lc_primary.h5',
        help='ZTF r-band primary HDF5 file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='data',
        help='Output directory'
    )

    args = parser.parse_args()

    # Expand paths
    ztf_g_path = os.path.expanduser(args.ztf_g)
    ztf_r_path = os.path.expanduser(args.ztf_r)

    print("="*70)
    print("Extract ZTF Lightcurves for Gaia Sources")
    print("="*70)

    # Step 1: Load Gaia SDSS names
    print("\nStep 1: Loading Gaia source list...")
    gaia = Table.read(args.gaia_file)
    gaia_sdss_names = set(str(name).strip() for name in gaia['sdss_name'])
    print(f"  Gaia sources: {len(gaia_sdss_names):,}")

    # Step 2: Extract ZTF g-band
    print("\nStep 2: Extracting ZTF g-band...")
    g_output = os.path.join(args.output_dir, 'ztf_g_for_gaia_primary.fits')
    extract_ztf_band(ztf_g_path, gaia_sdss_names, 'g', g_output)

    # Step 3: Extract ZTF r-band
    print("\nStep 3: Extracting ZTF r-band...")
    r_output = os.path.join(args.output_dir, 'ztf_r_for_gaia_primary.fits')
    extract_ztf_band(ztf_r_path, gaia_sdss_names, 'r', r_output)

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  {g_output}")
    print(f"  {r_output}")
    print("="*70)


if __name__ == '__main__':
    main()
