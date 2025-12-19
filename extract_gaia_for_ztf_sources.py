#!/usr/bin/env python
"""
Extract Gaia G-band lightcurves for SDSS QSOs that have ZTF primary lightcurves.

This script:
1. Loads ZTF primary HDF5 catalog to get SDSS names (731k sources)
2. Loads full Gaia source catalog (489k sources) to match sdss_name → source_id
3. Loads full Gaia epoch photometry (223k sources with lightcurves)
4. Extracts matching Gaia G-band lightcurves and saves to FITS

Usage:
    python extract_gaia_for_ztf_sources.py [--output OUTPUT.fits]
"""

import numpy as np
from astropy.table import Table
import argparse
import os


def main():
    parser = argparse.ArgumentParser(
        description='Extract Gaia G-band lightcurves for ZTF primary sources'
    )
    parser.add_argument(
        '--ztf-catalog',
        type=str,
        default='data/ztf_primary_sources_catalog.fits',
        help='ZTF primary sources catalog with sdss_name column'
    )
    parser.add_argument(
        '--gaia-source',
        type=str,
        default='~/data/gaia/sdssdr16q_gaia_source.fits',
        help='Full Gaia source catalog (sdss_name → source_id)'
    )
    parser.add_argument(
        '--gaia-epoch',
        type=str,
        default='~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits',
        help='Full Gaia epoch photometry (source_id → lightcurves)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='data/gaia_g_for_ztf_primary.fits',
        help='Output FITS file with Gaia G-band lightcurves'
    )

    args = parser.parse_args()

    # Expand paths
    gaia_source_path = os.path.expanduser(args.gaia_source)
    gaia_epoch_path = os.path.expanduser(args.gaia_epoch)

    print("="*70)
    print("Extract Gaia G-band for ZTF Primary Sources")
    print("="*70)

    # Step 1: Load ZTF primary SDSS names
    print("\nStep 1: Loading ZTF primary catalog...")
    ztf_cat = Table.read(args.ztf_catalog)
    ztf_sdss_names = set(str(name).strip() for name in ztf_cat['sdss_name'])
    print(f"  ZTF sources: {len(ztf_sdss_names):,}")

    # Step 2: Load Gaia source catalog to map sdss_name → source_id
    print("\nStep 2: Loading Gaia source catalog...")
    print(f"  Reading: {gaia_source_path}")
    gaia_src = Table.read(gaia_source_path, memmap=True)
    print(f"  Gaia sources: {len(gaia_src):,}")

    # Build mapping: sdss_name → (source_id, ra, dec, z)
    print("  Building sdss_name → source_id mapping...")
    sdss_to_gaia = {}
    for row in gaia_src:
        sdss_name = str(row['sdss_name']).strip()
        if sdss_name in ztf_sdss_names:
            source_id = row['source_id']
            ra = row['q_ra']
            dec = row['q_dec']
            z = row['z_qn'] if 'z_qn' in gaia_src.colnames else np.nan
            sdss_to_gaia[sdss_name] = (source_id, ra, dec, z)

    n_matched = len(sdss_to_gaia)
    print(f"  Matched: {n_matched:,} sources ({100*n_matched/len(ztf_sdss_names):.1f}% of ZTF)")

    if n_matched == 0:
        print("ERROR: No matching sources found!")
        return 1

    # Get set of matched source_ids
    matched_source_ids = set(sid for sid, _, _, _ in sdss_to_gaia.values())

    # Step 3: Load Gaia epoch photometry and filter to matched sources
    print("\nStep 3: Loading Gaia epoch photometry...")
    print(f"  Reading: {gaia_epoch_path}")
    gaia_epoch = Table.read(gaia_epoch_path, memmap=True)
    print(f"  Gaia lightcurves: {len(gaia_epoch):,}")

    # Filter to matched source_ids
    print("  Filtering to matched sources...")
    mask = np.isin(gaia_epoch['source_id'], list(matched_source_ids))
    gaia_matched = gaia_epoch[mask]
    print(f"  Kept: {len(gaia_matched):,} sources with lightcurves")

    # Step 4: Create reverse mapping: source_id → sdss_name
    print("\nStep 4: Creating output table...")
    source_id_to_sdss = {sid: name for name, (sid, _, _, _) in sdss_to_gaia.items()}

    # Build output arrays
    output_sdss_names = []
    output_source_ids = []
    output_ras = []
    output_decs = []
    output_zs = []

    for row in gaia_matched:
        sid = row['source_id']
        if sid in source_id_to_sdss:
            sdss_name = source_id_to_sdss[sid]
            _, ra, dec, z = sdss_to_gaia[sdss_name]

            output_sdss_names.append(sdss_name)
            output_source_ids.append(sid)
            output_ras.append(ra)
            output_decs.append(dec)
            output_zs.append(z)

    # Create output table with all Gaia columns plus sdss_name
    output = Table()
    output['sdss_name'] = output_sdss_names
    output['source_id'] = output_source_ids
    output['ra'] = output_ras
    output['dec'] = output_decs
    output['z'] = output_zs

    # Copy all lightcurve columns from gaia_matched
    for col in gaia_matched.colnames:
        if col != 'source_id':  # Already added
            output[col] = gaia_matched[col]

    # Step 5: Write output
    print(f"\nStep 5: Writing output...")
    output.write(args.output, format='fits', overwrite=True)

    file_size_mb = os.path.getsize(args.output) / (1024**2)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"ZTF sources: {len(ztf_sdss_names):,}")
    print(f"Gaia matches: {n_matched:,} ({100*n_matched/len(ztf_sdss_names):.1f}%)")
    print(f"With lightcurves: {len(output):,} ({100*len(output)/n_matched:.1f}% of matches)")
    print(f"\nOutput: {args.output}")
    print(f"Size: {file_size_mb:.1f} MB")
    print(f"Columns: {', '.join(output.colnames[:10])} ...")
    print("="*70)

    return 0


if __name__ == '__main__':
    exit(main())
