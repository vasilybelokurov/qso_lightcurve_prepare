#!/usr/bin/env python
"""
Merge ZTF chunk files with dtype correction for variable-length arrays.

Fixes the issue where empty arrays in airmass columns have dtype '<U1'
instead of 'float32', which prevents FITS file creation.

Adapted from qso_gp_mock/merge_ztf_chunks_fix_dtypes.py
"""
import os
import sys
import glob
import numpy as np
from astropy.table import Table, vstack
import argparse


def fix_column_dtypes(table, column_patterns):
    """
    Fix dtypes for columns matching patterns.

    Ensures all arrays in variable-length columns have consistent float32 dtype.
    """
    for col_pattern in column_patterns:
        matching_cols = [c for c in table.colnames if col_pattern in c]

        for col in matching_cols:
            # Convert each array in the column to float32
            fixed_data = []
            for val in table[col]:
                arr = np.array(val, dtype=np.float32)
                fixed_data.append(arr)
            table[col] = fixed_data

    return table


def merge_chunks(chunks_dir, output_file, column_patterns=['airmass', 'sharp', 'chi', 'limitmag']):
    """
    Merge chunk files with dtype correction.
    """
    chunk_files = sorted(glob.glob(os.path.join(chunks_dir, 'chunk_*.fits')))

    if len(chunk_files) == 0:
        print(f"ERROR: No chunk files found in {chunks_dir}")
        return 1

    print(f"Found {len(chunk_files)} chunk files")
    print(f"Fixing dtypes for columns: {column_patterns}")
    print()

    tables = []
    for i, chunk_file in enumerate(chunk_files):
        print(f"[{i+1}/{len(chunk_files)}] Processing {os.path.basename(chunk_file)}...")
        chunk = Table.read(chunk_file)

        # Show first chunk info
        if i == 0:
            print(f"  Columns: {chunk.colnames}")
            print(f"  Rows in first chunk: {len(chunk)}")
            print()

        # Fix dtypes
        chunk = fix_column_dtypes(chunk, column_patterns)

        tables.append(chunk)

    print(f"\nMerging {len(tables)} tables...")
    merged = vstack(tables)

    # Get file size estimate
    print(f"Total rows: {len(merged)}")
    print(f"Total columns: {len(merged.colnames)}")

    print(f"\nWriting to {output_file}...")
    merged.write(output_file, format='fits', overwrite=True)

    # Get actual file size
    file_size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"✓ Merged {len(merged)} sources to {output_file}")
    print(f"  File size: {file_size_mb:.1f} MB")

    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Merge ZTF chunks with dtype correction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--chunks-dir',
        type=str,
        default='ztf_chunks',
        help='Directory with chunk files'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='data/ztf_lightcurves_merged.fits',
        help='Output merged FITS file'
    )

    args = parser.parse_args()

    return merge_chunks(args.chunks_dir, args.output)


if __name__ == '__main__':
    sys.exit(main())
