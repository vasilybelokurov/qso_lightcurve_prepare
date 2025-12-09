#!/usr/bin/env python
"""
Apply Leave-One-Out (LOO) outlier cleaning to ZTF lightcurves.

Based on qso_gp_mock/ztf_10k_lightcurve_quality.ipynb cell 13.

Algorithm:
1. Fit DRW/OU GP to cleaned lightcurve
2. For each epoch, refit with that epoch removed (LOO)
3. Find epoch with largest |ΔlogD|
4. Remove if: ΔlogD < 0 (D decreases) AND |ΔlogD| >= 0.5
5. Repeat up to MAX_PASSES times

Usage:
    # Test on 100 sources
    python apply_loo_cleaning.py --test --n-test 100 --n-jobs 4

    # Full run
    python apply_loo_cleaning.py --n-jobs 8
"""

import os
import sys
import numpy as np
from astropy.table import Table
from joblib import Parallel, delayed
import argparse

from fit_drw import fit_drw_unconstrained


# LOO parameters (from notebook)
MIN_POINTS = 10
ONLY_IF_D_DECREASES = True
MIN_ABS_DELTA_LOGD = 0.5
MAX_FAILED_FRACTION = 0.10
MAX_LOO_PASSES = 5


def prepare_ztf_lightcurve(ztf_row, z, band='g', min_points=MIN_POINTS):
    """
    Extract and prepare ZTF lightcurve for fitting.

    Parameters
    ----------
    ztf_row : astropy.table.Row
        Row from ZTF lightcurve table
    z : float
        Redshift for rest-frame time conversion
    band : str
        Band to extract: 'g' or 'r'
    min_points : int
        Minimum number of valid epochs required

    Returns
    -------
    dict or None
        Dictionary with 'times_rest', 'mags', 'mag_errs' or None if unusable
    """
    band = band.lower()

    mjd_col = f'mjd_{band}'
    mag_col = f'mag_{band}'
    magerr_col = f'magerr_{band}'

    try:
        times = np.asarray(ztf_row[mjd_col], dtype=float)
        mags = np.asarray(ztf_row[mag_col], dtype=float)
        mag_errs = np.asarray(ztf_row[magerr_col], dtype=float)
    except (KeyError, ValueError):
        return None

    # Filter invalid values
    mask = ~(np.isnan(times) | np.isnan(mags) | np.isnan(mag_errs))
    mask &= (mag_errs > 0)

    if mask.sum() < min_points:
        return None

    times = times[mask]
    mags = mags[mask]
    mag_errs = mag_errs[mask]

    # Sort by time (celerite2 requires sorted times)
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    mags = mags[sort_idx]
    mag_errs = mag_errs[sort_idx]

    # Convert to rest-frame time
    times_rest = times / (1.0 + z)

    return {
        'times_rest': times_rest,
        'mags': mags,
        'mag_errs': mag_errs
    }


def _ztf_loo_process_one(i, ztf_table, catalog_table, band, min_points, max_passes):
    """
    Process a single source for LOO outlier detection.

    Returns
    -------
    tuple or None
        (result_dict, removed_orig, n_full_fail, n_loo_fail, n_post_fail, source_id)
        or None if source is skipped
    """
    ztf_row = ztf_table[i]
    cat_row = catalog_table[i]

    z = float(cat_row["z"])
    source_id = int(ztf_row["source_id"])

    lc = prepare_ztf_lightcurve(ztf_row, z, band=band, min_points=min_points)
    if lc is None:
        return None

    times_rest = np.asarray(lc["times_rest"], dtype=float)
    mags = np.asarray(lc["mags"], dtype=float)
    mag_errs = np.asarray(lc["mag_errs"], dtype=float)

    if len(times_rest) < min_points:
        return None

    mags_centered = mags - np.nanmedian(mags)

    times_cur = times_rest.copy()
    mags_cur = mags_centered.copy()
    errs_cur = mag_errs.copy()
    orig_indices = np.arange(len(times_rest), dtype=int)

    first_full_fit = None
    last_post_fit = None
    n_removed = 0
    n_passes_done = 0
    removed_orig = []

    n_full_fail = 0
    n_loo_fail = 0
    n_post_fail = 0

    while True:
        N_cur = len(times_cur)
        if N_cur < min_points:
            break

        # Full fit
        full_fit = fit_drw_unconstrained(times_cur, mags_cur, errs_cur)
        if full_fit is None or (not full_fit.get("success", True)):
            n_full_fail += 1
            last_post_fit = None
            break

        if first_full_fit is None:
            first_full_fit = full_fit

        log_sigma_full_cur = full_fit["log_sigma"]
        log_tau_full_cur = full_fit["log_tau"]
        logD_full_cur = 2.0 * log_sigma_full_cur - log_tau_full_cur

        # LOO fits
        delta_logD = np.full(N_cur, np.nan, dtype=float)
        n_failed_subfits = 0

        for j in range(N_cur):
            mask = np.ones(N_cur, dtype=bool)
            mask[j] = False
            fit_j = fit_drw_unconstrained(times_cur[mask], mags_cur[mask], errs_cur[mask])
            if fit_j is None or (not fit_j.get("success", True)):
                n_failed_subfits += 1
                continue

            log_sigma_j = fit_j["log_sigma"]
            log_tau_j = fit_j["log_tau"]
            logD_j = 2.0 * log_sigma_j - log_tau_j
            delta_logD[j] = logD_j - logD_full_cur

        frac_failed = n_failed_subfits / float(N_cur)
        if frac_failed > MAX_FAILED_FRACTION or np.all(~np.isfinite(delta_logD)):
            n_loo_fail += 1
            last_post_fit = None
            break

        idx_max_D_cur = int(np.nanargmax(np.abs(delta_logD)))
        delta_best = float(delta_logD[idx_max_D_cur])

        remove_point = True
        if ONLY_IF_D_DECREASES and not (delta_best < 0.0):
            remove_point = False
        if np.abs(delta_best) < MIN_ABS_DELTA_LOGD:
            remove_point = False

        if not remove_point:
            last_post_fit = full_fit
            break

        # Refit after removing most influential point
        mask_keep = np.ones(N_cur, dtype=bool)
        mask_keep[idx_max_D_cur] = False

        post_fit = fit_drw_unconstrained(times_cur[mask_keep], mags_cur[mask_keep], errs_cur[mask_keep])
        if post_fit is None or (not post_fit.get("success", True)):
            n_post_fail += 1
            last_post_fit = None
            break

        removed_orig.append(int(orig_indices[idx_max_D_cur]))

        times_cur = times_cur[mask_keep]
        mags_cur = mags_cur[mask_keep]
        errs_cur = errs_cur[mask_keep]
        orig_indices = orig_indices[mask_keep]

        last_post_fit = post_fit
        n_removed += 1
        n_passes_done += 1

        if n_passes_done >= max_passes:
            break

    if first_full_fit is None:
        return None

    # If LOO never accepted any removal, use full fit as post fit
    if last_post_fit is None:
        last_post_fit = first_full_fit

    log_sigma_full = first_full_fit["log_sigma"]
    log_tau_full = first_full_fit["log_tau"]
    logD_full = 2.0 * log_sigma_full - log_tau_full

    log_sigma_post = last_post_fit["log_sigma"]
    log_tau_post = last_post_fit["log_tau"]
    logD_post = 2.0 * log_sigma_post - log_tau_post

    result_dict = dict(
        source_id=source_id,
        z=z,
        band=band,
        log_sigma_full=log_sigma_full,
        log_tau_full=log_tau_full,
        logD_full=logD_full,
        log_sigma_post=log_sigma_post,
        log_tau_post=log_tau_post,
        logD_post=logD_post,
        n_epochs=len(times_rest),
        n_removed=n_removed,
        n_passes_done=n_passes_done,
    )

    return (result_dict, removed_orig, n_full_fail, n_loo_fail, n_post_fail, source_id)


def run_ztf_loo_for_band(ztf_table, catalog_table, band="g", min_points=MIN_POINTS,
                          max_passes=MAX_LOO_PASSES, n_jobs=-1):
    """
    Run LOO outlier detection for one band using parallel processing.

    Parameters
    ----------
    ztf_table : astropy.table.Table
        ZTF lightcurve table
    catalog_table : astropy.table.Table
        Catalog table with z column
    band : str
        Band to process ('g' or 'r')
    min_points : int
        Minimum epochs required
    max_passes : int
        Maximum LOO removal passes
    n_jobs : int
        Number of parallel jobs (-1 = all cores)

    Returns
    -------
    results_table : astropy.table.Table
        Table with LOO results
    removed_indices : dict
        Dictionary mapping source_id -> array of removed epoch indices
    """
    if n_jobs is None or n_jobs == -1:
        n_jobs = -1

    n_obj = len(ztf_table)
    indices = np.arange(n_obj)

    results = []
    removed_indices = {}
    total_full_fail = 0
    total_loo_fail = 0
    total_post_fail = 0

    def wrapped(i):
        return _ztf_loo_process_one(i, ztf_table, catalog_table, band, min_points, max_passes)

    print(f"Running ZTF LOO for band {band} with joblib on {n_obj} objects, n_jobs={n_jobs}")

    parallel = Parallel(n_jobs=n_jobs, prefer="processes", return_as="generator")
    gen = parallel(delayed(wrapped)(i) for i in indices)

    count = 0
    log_interval = 100  # Log every 100 sources
    for out in gen:
        if out is None:
            continue
        (res_dict, removed_orig, n_full_fail, n_loo_fail, n_post_fail, source_id) = out

        results.append(res_dict)
        removed_indices[source_id] = np.array(removed_orig, dtype=int)
        total_full_fail += n_full_fail
        total_loo_fail += n_loo_fail
        total_post_fail += n_post_fail

        count += 1
        if count % log_interval == 0:
            print(f"  Processed {count}/{n_obj} sources ({100*count/n_obj:.1f}%)", flush=True)

    results_table = Table(rows=results)

    print(f"Finished ZTF LOO for band {band}")
    print(f"  Successful fits: {len(results)}/{n_obj}")
    print(f"  Full-fit failures: {total_full_fail}")
    print(f"  LOO failures: {total_loo_fail}")
    print(f"  Post-fit failures: {total_post_fail}")

    return results_table, removed_indices


def create_loo_cleaned_lightcurves(ztf_table, removed_dict_g, removed_dict_r, bands=['g', 'r']):
    """
    Create cleaned ZTF lightcurve table with LOO-identified outliers removed.

    Parameters
    ----------
    ztf_table : astropy.table.Table
        Original ZTF lightcurve table
    removed_dict_g : dict
        Dictionary mapping source_id -> removed indices for g-band
    removed_dict_r : dict
        Dictionary mapping source_id -> removed indices for r-band
    bands : list
        Bands to process

    Returns
    -------
    astropy.table.Table
        Cleaned ZTF table with outliers removed
    """
    print("\nCreating LOO-cleaned lightcurve table...")

    ztf_cleaned = Table()
    ztf_cleaned['source_id'] = ztf_table['source_id']
    ztf_cleaned['ra'] = ztf_table['ra']
    ztf_cleaned['dec'] = ztf_table['dec']

    removed_dicts = {'g': removed_dict_g, 'r': removed_dict_r}

    for band in bands:
        print(f"  Processing {band}-band...")
        removed_dict = removed_dicts[band]

        mjd_clean = []
        mag_clean = []
        magerr_clean = []
        catflags_clean = []
        sharp_clean = []
        chi_clean = []
        limitmag_clean = []
        airmass_clean = []

        n_removed_total = 0

        for i in range(len(ztf_table)):
            source_id = int(ztf_table['source_id'][i])

            # Get removed indices for this source
            if source_id in removed_dict:
                removed_idx = removed_dict[source_id]
                n_removed_total += len(removed_idx)
            else:
                removed_idx = np.array([], dtype=int)

            # Create keep mask
            n_epochs = len(ztf_table[f'mjd_{band}'][i])
            keep_mask = np.ones(n_epochs, dtype=bool)
            if len(removed_idx) > 0:
                keep_mask[removed_idx] = False

            # Apply mask to all columns
            mjd_clean.append(np.array(ztf_table[f'mjd_{band}'][i][keep_mask], dtype=float))
            mag_clean.append(np.array(ztf_table[f'mag_{band}'][i][keep_mask], dtype=np.float32))
            magerr_clean.append(np.array(ztf_table[f'magerr_{band}'][i][keep_mask], dtype=np.float32))
            catflags_clean.append(np.array(ztf_table[f'catflags_{band}'][i][keep_mask], dtype=np.int32))
            sharp_clean.append(np.array(ztf_table[f'sharp_{band}'][i][keep_mask], dtype=np.float32))
            chi_clean.append(np.array(ztf_table[f'chi_{band}'][i][keep_mask], dtype=np.float32))
            limitmag_clean.append(np.array(ztf_table[f'limitmag_{band}'][i][keep_mask], dtype=np.float32))
            airmass_clean.append(np.array(ztf_table[f'airmass_{band}'][i][keep_mask], dtype=np.float32))

        # Add to cleaned table
        ztf_cleaned[f'mjd_{band}'] = mjd_clean
        ztf_cleaned[f'mag_{band}'] = mag_clean
        ztf_cleaned[f'magerr_{band}'] = magerr_clean
        ztf_cleaned[f'catflags_{band}'] = catflags_clean
        ztf_cleaned[f'sharp_{band}'] = sharp_clean
        ztf_cleaned[f'chi_{band}'] = chi_clean
        ztf_cleaned[f'limitmag_{band}'] = limitmag_clean
        ztf_cleaned[f'airmass_{band}'] = airmass_clean

        print(f"    {band}-band: removed {n_removed_total} epochs total")

    return ztf_cleaned


def main():
    parser = argparse.ArgumentParser(
        description="Apply LOO outlier cleaning to ZTF lightcurves",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-ztf', type=str, default='data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits',
                       help='Input cleaned ZTF lightcurves')
    parser.add_argument('--input-catalog', type=str, default='data/gaia_ztf_qso_sample_catalog.fits',
                       help='Input catalog (for z)')
    parser.add_argument('--output-results', type=str, default=None,
                       help='Output LOO results table (auto-generated if not provided)')
    parser.add_argument('--output-cleaned', type=str, default=None,
                       help='Output LOO-cleaned lightcurves (auto-generated if not provided)')
    parser.add_argument('--band', type=str, required=True, choices=['g', 'r'],
                       help='Band to process (g or r) - REQUIRED')
    parser.add_argument('--n-jobs', type=int, default=-1, help='Number of parallel jobs (-1=all)')
    parser.add_argument('--test', action='store_true', help='Test mode (use subset)')
    parser.add_argument('--n-test', type=int, default=100, help='Number of sources for test mode')

    args = parser.parse_args()

    # Single band processing
    band = args.band.lower()
    bands = [band]

    # Auto-generate output filenames if not provided
    if args.output_results is None:
        args.output_results = f'data/ztf_loo_results_{band}.fits'
    if args.output_cleaned is None:
        args.output_cleaned = f'data/gaia_ztf_qso_sample_ztf_lc_{band}_clean_loo.fits'

    print("="*70)
    print("ZTF LOO OUTLIER CLEANING (SINGLE BAND)")
    print("="*70)
    print(f"Input ZTF:     {args.input_ztf}")
    print(f"Input catalog: {args.input_catalog}")
    print(f"Output results: {args.output_results}")
    print(f"Output cleaned: {args.output_cleaned}")
    print(f"Band:          {band}")
    print(f"Parallel jobs: {args.n_jobs}")
    if args.test:
        print(f"TEST MODE:     Using first {args.n_test} sources")
    print("="*70)

    # Load data
    print("\nLoading data...")
    ztf = Table.read(args.input_ztf)
    catalog = Table.read(args.input_catalog)
    print(f"  ZTF: {len(ztf)} sources")
    print(f"  Catalog: {len(catalog)} sources")

    # Verify alignment
    assert np.all(ztf['source_id'] == catalog['source_id']), "ZTF-Catalog alignment failed!"
    print("  ✓ ZTF-Catalog aligned")

    # Test mode subset
    if args.test:
        ztf = ztf[:args.n_test]
        catalog = catalog[:args.n_test]
        print(f"\nTest mode: processing {len(ztf)} sources")

    # Run LOO for single band
    print(f"\n{'='*70}")
    print(f"Processing {band}-band")
    print("="*70)

    results_table, removed_dict = run_ztf_loo_for_band(
        ztf, catalog, band=band, n_jobs=args.n_jobs
    )

    # Statistics
    n_removed_arr = [len(removed_dict[sid]) for sid in removed_dict.keys()]
    print(f"\n{band}-band LOO statistics:")
    print(f"  Mean removed: {np.mean(n_removed_arr):.2f}")
    print(f"  Median removed: {np.median(n_removed_arr):.0f}")
    print(f"  Max removed: {np.max(n_removed_arr)}")
    print(f"  No removal: {sum(1 for n in n_removed_arr if n == 0)}/{len(n_removed_arr)}")

    # Save results
    print(f"\n{'='*70}")
    print("SAVING RESULTS")
    print("="*70)

    results_table.write(args.output_results, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(args.output_results) / (1024 * 1024)
    print(f"✓ Saved LOO results: {args.output_results}")
    print(f"  {len(results_table)} rows ({file_size_mb:.1f} MB)")

    # Create cleaned lightcurves for single band
    ztf_cleaned = create_loo_cleaned_lightcurves(
        ztf, removed_dict if band == 'g' else {}, removed_dict if band == 'r' else {}, bands
    )

    ztf_cleaned.write(args.output_cleaned, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(args.output_cleaned) / (1024 * 1024)
    print(f"\n✓ Saved LOO-cleaned lightcurves: {args.output_cleaned}")
    print(f"  {len(ztf_cleaned)} sources ({file_size_mb:.1f} MB)")

    print("\n" + "="*70)
    print("LOO CLEANING COMPLETE")
    print("="*70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
