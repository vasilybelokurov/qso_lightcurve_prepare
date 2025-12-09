#!/usr/bin/env python3
"""
Create test catalog with first 2000 sources for ZTF download testing.
"""

from astropy.table import Table
import os

def create_test_catalog():
    """Extract first 2000 sources from full target catalog."""

    print("Creating test catalog (first 2000 sources)...")

    # Load full catalog
    input_file = 'data/sdssdr16q_gaia_ztf_targets.fits'
    catalog = Table.read(input_file)

    print(f"  Full catalog: {len(catalog)} sources")

    # Extract first 2000
    test_catalog = catalog[:2000]

    # Save
    output_file = 'data/test_2k_ztf_targets.fits'
    test_catalog.write(output_file, overwrite=True)

    file_size_kb = os.path.getsize(output_file) / 1024

    print(f"  Test catalog: {len(test_catalog)} sources")
    print(f"  Output: {output_file}")
    print(f"  Size: {file_size_kb:.1f} KB")
    print()
    print("Ready for testing!")

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(create_test_catalog())
