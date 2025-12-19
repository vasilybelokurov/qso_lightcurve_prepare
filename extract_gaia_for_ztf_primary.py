#!/usr/bin/env python
"""
Extract Gaia G-band lightcurves for SDSS QSOs that have ZTF primary lightcurves.

This script:
1. Loads ZTF primary HDF5 files (g and r bands) to get the list of SDSS names
2. Loads Gaia catalog and LOO-cleaned lightcurves
3. Extracts Gaia G-band lightcurves for matching SDSS names
4. Saves to HDF5 file in the same format as ZTF primary files

Output: gaia_g_band_lc_for_ztf_primary.h5

Usage:
    python extract_gaia_for_ztf_primary.py [--output OUTPUT.h5]
"""

import h5py
import numpy as np
from astropy.table import Table
import argparse
import os
from tqdm import tqdm


def load_ztf_sdss_names(ztf_g_path: str, ztf_r_path: str) -> set:
    """
    Load all unique SDSS names from ZTF primary g and r band files.

    Parameters
    ----------
    ztf_g_path : str
        Path to ZTF g-band primary HDF5 file
    ztf_r_path : str
        Path to ZTF r-band primary HDF5 file

    Returns
    -------
    set
        Set of SDSS names (strings)
    """
    print("Loading ZTF primary SDSS names...")

    sdss_names = set()

    # Load g-band
    with h5py.File(ztf_g_path, 'r') as f:
        for name in f['sdss_name'][:]:
            if isinstance(name, bytes):
                sdss_names.add(name.decode().strip())
            else:
                sdss_names.add(str(name).strip())

    print(f"  ZTF g-band: {len(sdss_names)} unique SDSS names")

    # Load r-band
    with h5py.File(ztf_r_path, 'r') as f:
        for name in f['sdss_name'][:]:
            if isinstance(name, bytes):
                sdss_names.add(name.decode().strip())
            else:
                sdss_names.add(str(name).strip())

    print(f"  ZTF g+r combined: {len(sdss_names)} unique SDSS names")

    return sdss_names


def extract_gaia_g_band(
    gaia_catalog_path: str,
    gaia_lc_path: str,
    ztf_sdss_names: set,
    output_path: str
):
    """
    Extract Gaia G-band lightcurves for sources matching ZTF primary list.

    Parameters
    ----------
    gaia_catalog_path : str
        Path to Gaia catalog FITS file
    gaia_lc_path : str
        Path to Gaia G-band LOO-cleaned lightcurves FITS file
    ztf_sdss_names : set
        Set of SDSS names from ZTF primary files
    output_path : str
        Output HDF5 file path
    """
    print("\nLoading Gaia data...")

    # Load Gaia catalog
    cat = Table.read(gaia_catalog_path)
    print(f"  Gaia catalog: {len(cat)} sources")

    # Load Gaia G-band LOO-cleaned lightcurves
    gaia_lc = Table.read(gaia_lc_path)
    print(f"  Gaia G-band LOO cleaned: {len(gaia_lc)} sources")

    # Match by SDSS name
    print("\nMatching sources...")

    matched_indices = []
    matched_sdss_names = []

    for i, name in enumerate(tqdm(cat['sdss_name'], desc="  Matching")):
        name_str = str(name).strip()
        if name_str in ztf_sdss_names:
            matched_indices.append(i)
            matched_sdss_names.append(name_str)

    n_matched = len(matched_indices)
    print(f"  Matched: {n_matched} sources ({100*n_matched/len(ztf_sdss_names):.1f}% of ZTF)")

    if n_matched == 0:
        print("ERROR: No matching sources found!")
        return

    # Extract matched data
    print("\nExtracting Gaia G-band lightcurves...")

    matched_cat = cat[matched_indices]
    matched_lc = gaia_lc[matched_indices]

    # Prepare output data
    output_data = {
        'sdss_name': [],
        'source_id': [],
        'ra': [],
        'dec': [],
        'z': [],
        'mjd': [],
        'mag': [],
        'magerr': [],
        'flux': [],
        'flux_error': [],
        'n_epochs': [],
    }

    for i in tqdm(range(n_matched), desc="  Processing"):
        # Get SDSS info from catalog
        sdss_name = matched_sdss_names[i]
        source_id = matched_lc['source_id'][i]
        ra = matched_cat['q_ra'][i]
        dec = matched_cat['q_dec'][i]
        z = matched_cat['z_qn'][i] if 'z_qn' in matched_cat.colnames else np.nan

        # Get Gaia G-band lightcurve (already LOO-cleaned)
        mjd = matched_lc['g_transit_time'][i]
        mag = matched_lc['g_transit_mag'][i]
        flux = matched_lc['g_transit_flux'][i]
        flux_error = matched_lc['g_transit_flux_error'][i]

        # Remove NaN values (already filtered by LOO cleaning, but be safe)
        valid = ~(np.isnan(mjd) | np.isnan(mag) | np.isnan(flux))
        mjd = mjd[valid]
        mag = mag[valid]
        flux = flux[valid]
        flux_error = flux_error[valid]

        # Compute magnitude error from flux error
        # magerr = 2.5/ln(10) * flux_err/flux = 1.0857 * flux_err/flux
        with np.errstate(divide='ignore', invalid='ignore'):
            magerr = 1.0857 * flux_error / np.abs(flux)
            magerr = np.where(np.isfinite(magerr), magerr, np.nan)

        n_epochs = len(mjd)

        # Store
        output_data['sdss_name'].append(sdss_name)
        output_data['source_id'].append(source_id)
        output_data['ra'].append(ra)
        output_data['dec'].append(dec)
        output_data['z'].append(z)
        output_data['mjd'].append(mjd)
        output_data['mag'].append(mag)
        output_data['magerr'].append(magerr)
        output_data['flux'].append(flux)
        output_data['flux_error'].append(flux_error)
        output_data['n_epochs'].append(n_epochs)

    # Write to HDF5
    print(f"\nWriting to {output_path}...")

    with h5py.File(output_path, 'w') as f:
        # Metadata
        f.attrs['band'] = 'G'
        f.attrs['n_sources'] = n_matched
        f.attrs['description'] = 'Gaia G-band lightcurves for SDSS QSOs with ZTF primary lightcurves'
        f.attrs['source'] = 'Gaia DR3 epoch photometry (LOO-cleaned)'
        f.attrs['created_by'] = 'extract_gaia_for_ztf_primary.py'

        # Fixed-length columns
        f.create_dataset('sdss_name',
                        data=np.array(output_data['sdss_name'], dtype='S19'),
                        compression='gzip', compression_opts=9)
        f.create_dataset('source_id',
                        data=np.array(output_data['source_id'], dtype=np.int64),
                        compression='gzip', compression_opts=9)
        f.create_dataset('ra',
                        data=np.array(output_data['ra'], dtype=np.float64),
                        compression='gzip', compression_opts=9)
        f.create_dataset('dec',
                        data=np.array(output_data['dec'], dtype=np.float64),
                        compression='gzip', compression_opts=9)
        f.create_dataset('z',
                        data=np.array(output_data['z'], dtype=np.float32),
                        compression='gzip', compression_opts=9)
        f.create_dataset('n_epochs',
                        data=np.array(output_data['n_epochs'], dtype=np.int32),
                        compression='gzip', compression_opts=9)

        # Variable-length arrays
        vlen_float64 = h5py.vlen_dtype(np.float64)
        vlen_float32 = h5py.vlen_dtype(np.float32)

        f.create_dataset('mjd', data=output_data['mjd'], dtype=vlen_float64)
        f.create_dataset('mag', data=output_data['mag'], dtype=vlen_float32)
        f.create_dataset('magerr', data=output_data['magerr'], dtype=vlen_float32)
        f.create_dataset('flux', data=output_data['flux'], dtype=vlen_float32)
        f.create_dataset('flux_error', data=output_data['flux_error'], dtype=vlen_float32)

    # Statistics
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Matched sources: {n_matched}")
    print(f"Total epochs: {sum(output_data['n_epochs']):,}")

    n_epochs_array = np.array(output_data['n_epochs'])
    print(f"\nEpochs per source:")
    print(f"  Median: {np.median(n_epochs_array):.0f}")
    print(f"  Mean: {np.mean(n_epochs_array):.1f}")
    print(f"  Range: {n_epochs_array.min()}-{n_epochs_array.max()}")

    file_size_mb = os.path.getsize(output_path) / (1024**2)
    print(f"\nOutput file size: {file_size_mb:.1f} MB")
    print(f"Output: {output_path}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Extract Gaia G-band lightcurves for SDSS QSOs with ZTF primary lightcurves'
    )
    parser.add_argument(
        '--ztf-g',
        type=str,
        default='~/data/ztf/qso/ztf_g_band_lc_primary.h5',
        help='Path to ZTF g-band primary HDF5 file'
    )
    parser.add_argument(
        '--ztf-r',
        type=str,
        default='~/data/ztf/qso/ztf_r_band_lc_primary.h5',
        help='Path to ZTF r-band primary HDF5 file'
    )
    parser.add_argument(
        '--gaia-catalog',
        type=str,
        default='data/gaia_ztf_qso_sample_catalog.fits',
        help='Path to Gaia catalog FITS file'
    )
    parser.add_argument(
        '--gaia-lc',
        type=str,
        default='data/gaia_ztf_qso_sample_gaia_lc_g_clean_loo.fits',
        help='Path to Gaia G-band LOO-cleaned lightcurves FITS file'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='data/gaia_g_band_lc_for_ztf_primary.h5',
        help='Output HDF5 file path'
    )

    args = parser.parse_args()

    # Expand paths
    ztf_g_path = os.path.expanduser(args.ztf_g)
    ztf_r_path = os.path.expanduser(args.ztf_r)
    gaia_catalog_path = args.gaia_catalog
    gaia_lc_path = args.gaia_lc
    output_path = args.output

    print("="*60)
    print("Extract Gaia G-band for ZTF Primary Sources")
    print("="*60)
    print(f"\nInput files:")
    print(f"  ZTF g-band: {ztf_g_path}")
    print(f"  ZTF r-band: {ztf_r_path}")
    print(f"  Gaia catalog: {gaia_catalog_path}")
    print(f"  Gaia G-band LC: {gaia_lc_path}")
    print(f"\nOutput: {output_path}")
    print()

    # Load ZTF SDSS names
    ztf_sdss_names = load_ztf_sdss_names(ztf_g_path, ztf_r_path)

    # Extract Gaia G-band for matching sources
    extract_gaia_g_band(
        gaia_catalog_path,
        gaia_lc_path,
        ztf_sdss_names,
        output_path
    )

    print("\n✓ Complete")


if __name__ == '__main__':
    main()
