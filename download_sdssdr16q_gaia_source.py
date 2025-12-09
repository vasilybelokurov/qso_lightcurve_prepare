#!/usr/bin/env python3
"""
Download SDSS DR16 QSO × Gaia DR3 source cross-match catalog.

Cross-matches all SDSS DR16Q final QSOs with Gaia DR3 sources within 1 arcsec.
Saves all columns from both catalogs to FITS file.

Output: ~/data/gaia/sdssdr16q_gaia_source.fits
"""

import sqlutilpy as sqlutil
from astropy.table import Table
import os
import time


def get_column_list_with_prefix(schema, table, prefix, conflicts):
    """
    Generate column list with prefix for conflicting columns.

    Parameters
    ----------
    schema : str
        Database schema name
    table : str
        Table name
    prefix : str
        Prefix to add to conflicting columns
    conflicts : list
        List of column names that conflict

    Returns
    -------
    str
        SQL column list (e.g., "col1, col2 AS prefix_col3, ...")
    """
    # Get all columns for this table
    query = f"""
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema='{schema}' AND table_name='{table}'
    ORDER BY ordinal_position
    """

    result = sqlutil.get(query, asDict=True)
    columns = result['column_name']

    # Build SELECT clause
    select_items = []
    for col in columns:
        if col in conflicts:
            select_items.append(f"{prefix}.{col} AS {prefix}_{col}")
        else:
            select_items.append(f"{prefix}.{col}")

    return ",\n        ".join(select_items)


def download_sdssdr16q_gaia_crossmatch():
    """Download SDSS DR16Q × Gaia DR3 cross-match with all columns."""

    print("=" * 70)
    print("SDSS DR16Q × Gaia DR3 Source Cross-Match")
    print("=" * 70)
    print()

    # Define conflicting column names (present in both tables)
    conflicts = ['ra', 'dec', 'duplicated_source']

    print("Step 1: Generating SQL query with column renaming...")

    # Build column lists
    sdss_cols = get_column_list_with_prefix('sdssdr16qso', 'main', 'q', conflicts)
    gaia_cols = get_column_list_with_prefix('gaia_dr3', 'gaia_source', 'g', conflicts)

    # Build full query
    query = f"""
    SELECT
        -- SDSS DR16Q columns (conflicts prefixed with sdss_)
        {sdss_cols},

        -- Gaia DR3 columns (conflicts prefixed with gaia_)
        {gaia_cols}
    FROM sdssdr16qso.main AS q
    JOIN gaia_dr3.gaia_source AS g
      ON q3c_join(q.ra, q.dec, g.ra, g.dec, 1.0/3600.0)
    WHERE q.is_qso_final = 1
    """

    print(f"  Conflicting columns to be prefixed: {conflicts}")
    print()

    print("Step 2: Executing cross-match query (1 arcsec radius)...")
    print("  This may take several minutes...")

    t0 = time.time()
    data = sqlutil.get(query, asDict=True)
    elapsed = time.time() - t0

    n_matched = len(data[list(data.keys())[0]]) if data else 0

    print(f"  ✓ Cross-matched {n_matched} QSOs in {elapsed:.1f} seconds")
    print()

    # Convert to table
    result = Table(data)

    print(f"Step 3: Saving catalog...")
    output_file = os.path.expanduser('~/data/gaia/sdssdr16q_gaia_source.fits')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    result.write(output_file, overwrite=True)

    # Get file size
    file_size_mb = os.path.getsize(output_file) / (1024 * 1024)

    print(f"  ✓ Saved: {output_file}")
    print(f"  File size: {file_size_mb:.1f} MB")
    print(f"  Rows: {len(result)}")
    print(f"  Columns: {len(result.colnames)}")
    print()

    # Print summary statistics
    print("=" * 70)
    print("Summary Statistics")
    print("=" * 70)
    print(f"Total QSOs matched with Gaia: {len(result)}")

    # Check for source_id column
    if 'source_id' in result.colnames:
        print(f"Unique Gaia source_ids: {len(set(result['source_id']))}")

    # Check redshift distribution
    if 'z' in result.colnames:
        import numpy as np
        print(f"Redshift range: {np.min(result['z']):.3f} - {np.max(result['z']):.3f}")
        print(f"Median redshift: {np.median(result['z']):.3f}")

    print("=" * 70)
    print()
    print("Next step: Run download_sdssdr16q_gaia_epoch_photometry.py")

    return result


def main():
    result = download_sdssdr16q_gaia_crossmatch()
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
