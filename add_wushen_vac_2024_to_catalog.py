#!/usr/bin/env python3
"""
Augment an existing SDSS+Gaia catalog FITS with extra columns from WSDB.

Default use-case: take the ~220k matched catalog
`data/sdss_gaia_catalog_for_ztf_gaia_lc.fits` and add columns from:
    sdssdr16qso.wushen_vac_2024

This script uses sqlutilpy (WSDB) and a `local_join` upload to preserve the
exact row order of the input catalog.
"""

from __future__ import annotations

import argparse
import os
from typing import Iterable

import numpy as np
import sqlutilpy as sqlutil
from astropy.table import Table


def _qident(name: str) -> str:
    escaped = name.replace('"', '""')
    return f'"{escaped}"'


def _split_qualified_table(qualified: str) -> tuple[str, str]:
    if "." not in qualified:
        raise ValueError(f"Expected qualified table like 'schema.table', got: {qualified!r}")
    schema, table = qualified.split(".", 1)
    return schema, table


def fetch_table_columns(schema: str, table: str) -> list[str]:
    query = f"""
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema='{schema}' AND table_name='{table}'
    ORDER BY ordinal_position
    """
    result = sqlutil.get(query, asDict=True)
    cols = list(result["column_name"])
    if not cols:
        raise RuntimeError(f"No columns found for {schema}.{table} (check permissions/name).")
    return cols


def build_select_list(
    vac_columns: Iterable[str],
    *,
    vac_alias: str,
    input_columns: set[str],
    join_key_vac: str,
    prefix: str,
    prefix_all: bool,
) -> tuple[list[str], list[str]]:
    select_items: list[str] = []
    output_names: list[str] = []

    for col in vac_columns:
        if col == join_key_vac:
            continue
        out = f"{prefix}{col}" if (prefix_all or col in input_columns) else col
        select_items.append(f'{vac_alias}.{_qident(col)} AS {_qident(out)}')
        output_names.append(out)

    return select_items, output_names


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Add columns from sdssdr16qso.wushen_vac_2024 to an existing FITS catalog",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        default="data/sdss_gaia_catalog_for_ztf_gaia_lc.fits",
        help="Input catalog FITS (must contain the join key column)",
    )
    parser.add_argument(
        "--output",
        default="data/sdss_gaia_catalog_for_ztf_gaia_lc_with_wushen_vac_2024.fits",
        help="Output FITS path",
    )
    parser.add_argument(
        "--vac-table",
        default="sdssdr16qso.wushen_vac_2024",
        help="WSDB table to join in (schema.table)",
    )
    parser.add_argument(
        "--input-key",
        default="sdss_name",
        help="Join key column name in the input FITS",
    )
    parser.add_argument(
        "--vac-key",
        default="sdss_name",
        help="Join key column name in the VAC table",
    )
    parser.add_argument(
        "--prefix",
        default="wushen_",
        help="Prefix applied to conflicting column names (or all columns with --prefix-all)",
    )
    parser.add_argument(
        "--prefix-all",
        action="store_true",
        help="Prefix all VAC columns (not just conflicts) to keep provenance obvious",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the WSDB SQL and exit without downloading",
    )
    args = parser.parse_args()

    input_path = os.path.expanduser(args.input)
    output_path = os.path.expanduser(args.output)

    print("=" * 70)
    print("ADD WUSHEN VAC 2024 COLUMNS")
    print("=" * 70)
    print(f"Input:    {input_path}")
    print(f"Output:   {output_path}")
    print(f"VAC:      {args.vac_table}")
    print(f"Join:     input.{args.input_key} = vac.{args.vac_key}")
    print(f"Prefix:   {args.prefix} (prefix_all={args.prefix_all})")
    print("=" * 70)

    print("\nLoading input FITS (this can be large)...")
    cat = Table.read(input_path, memmap=True)
    n_rows = len(cat)
    print(f"  Rows: {n_rows:,}")
    print(f"  Columns: {len(cat.colnames)}")

    if args.input_key not in cat.colnames:
        raise SystemExit(f"ERROR: join key {args.input_key!r} not found in input columns.")

    schema, table = _split_qualified_table(args.vac_table)
    print(f"\nFetching VAC columns list for {schema}.{table}...")
    vac_columns = fetch_table_columns(schema, table)
    if args.vac_key not in vac_columns:
        raise SystemExit(
            f"ERROR: VAC join key {args.vac_key!r} not found in {schema}.{table} columns."
        )
    print(f"  VAC columns: {len(vac_columns)}")

    select_items, output_names = build_select_list(
        vac_columns,
        vac_alias="v",
        input_columns=set(cat.colnames),
        join_key_vac=args.vac_key,
        prefix=args.prefix,
        prefix_all=args.prefix_all,
    )
    if not select_items:
        raise SystemExit("ERROR: no VAC columns selected (unexpected).")

    # Upload an explicit index to preserve order deterministically.
    idx = np.arange(n_rows, dtype=np.int64)
    join_values_list = [str(x).strip() for x in cat[args.input_key]]
    max_len = max((len(s) for s in join_values_list), default=1)
    # sqlutilpy expects numpy arrays with a real dtype (object arrays break upload schema inference)
    join_values = np.array(join_values_list, dtype=f"U{max_len}")

    # DISTINCT ON (m.idx) protects against accidental 1-to-many joins in the VAC.
    sql = f"""
    SELECT DISTINCT ON (m.idx)
        m.idx,
        m.{_qident(args.input_key)} AS {_qident(args.input_key)},
        {", ".join(select_items)}
    FROM mytmptable AS m
    LEFT JOIN {_qident(schema)}.{_qident(table)} AS v
      ON v.{_qident(args.vac_key)} = m.{_qident(args.input_key)}
    ORDER BY m.idx
    """

    if args.dry_run:
        print("\n--- SQL (dry run) ---")
        print(sql.strip())
        return 0

    print("\nQuerying WSDB via local_join (requires WSDB credentials/network access)...")
    result = sqlutil.local_join(
        sql,
        "mytmptable",
        (idx, join_values),
        ("idx", args.input_key),
        asDict=True,
    )

    # Build a table from the result; keep idx for alignment checks.
    joined = Table(result)
    if len(joined) != n_rows:
        raise SystemExit(
            f"ERROR: join returned {len(joined):,} rows, expected {n_rows:,}. "
            "This suggests a non-unique join key in the VAC table."
        )
    if not np.all(np.asarray(joined["idx"]) == idx):
        raise SystemExit("ERROR: WSDB result row order mismatch; idx alignment failed.")
    if not np.all(np.asarray(joined[args.input_key]).astype(str) == join_values.astype(str)):
        raise SystemExit("ERROR: WSDB result join-key mismatch; alignment check failed.")

    print("\nAdding VAC columns to input table...")
    added = 0
    for name in output_names:
        if name in cat.colnames:
            raise SystemExit(
                f"ERROR: output column {name!r} already exists after conflict handling; "
                "try --prefix-all or choose a different --prefix."
            )
        cat[name] = joined[name]
        added += 1
    print(f"  Added columns: {added}")

    print(f"\nWriting output FITS: {output_path}")
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    cat.write(output_path, format="fits", overwrite=True)

    size_mb = os.path.getsize(output_path) / (1024**2)
    print(f"  ✓ Done: {size_mb:.1f} MB, {len(cat):,} rows, {len(cat.colnames)} cols")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
