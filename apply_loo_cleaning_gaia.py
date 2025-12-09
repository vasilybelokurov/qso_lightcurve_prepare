#!/usr/bin/env python
"""
Apply LOO (Leave-One-Out) outlier cleaning to Gaia lightcurves.

Fits DRW/OU GP models to Gaia G/BP/RP lightcurves and iteratively removes
epochs that significantly increase the damping parameter D when left out.

Based on notebook: qso_gp_mock/ztf_10k_lightcurve_quality.ipynb (cell 14)

Usage:
    python apply_loo_cleaning_gaia.py --band g --n-jobs 8
    python apply_loo_cleaning_gaia.py --band bp --n-jobs 8
    python apply_loo_cleaning_gaia.py --band rp --n-jobs 8
"""

import os
import sys
import numpy as np
from astropy.table import Table
import argparse
from joblib import Parallel, delayed

from fit_drw import fit_drw_unconstrained

# LOO parameters (from qso_gp_mock notebook)
MIN_POINTS = 10
MAX_LOO_PASSES = 5
ONLY_IF_D_DECREASES = True
MIN_ABS_DELTA_LOGD = 0.5
MAX_FAILED_FRACTION = 0.10

# Helpful constant for flux->mag error conversion
K_MAG = 2.5 / np.log(10.0)


def prepare_gaia_lightcurve(gaia_row, z, band='g', min_points=MIN_POINTS):
    """
    Extract and clean Gaia lightcurve for one band.

    Parameters
    ----------
    gaia_row : astropy.table.Row
        Row from Gaia lightcurve table
    z : float
        Redshift for rest-frame time conversion
    band : str
        Band to extract: 'g', 'bp', or 'rp'
    min_points : int
        Minimum number of valid epochs required

    Returns
    -------
    dict or None
        Dictionary with 'times_rest', 'mags', 'mag_errs' or None if unusable
    """
    band = band.lower()

    # Map band name to Gaia column names
    if band == 'g':
        time_col = 'g_transit_time'
        mag_col = 'g_transit_mag'
        flux_col = 'g_transit_flux'
        flux_err_col = 'g_transit_flux_error'
    elif band == 'bp':
        time_col = 'bp_obs_time'
        mag_col = 'bp_mag'
        flux_col = 'bp_flux'
        flux_err_col = 'bp_flux_error'
    elif band == 'rp':
        time_col = 'rp_obs_time'
        mag_col = 'rp_mag'
        flux_col = 'rp_flux'
        flux_err_col = 'rp_flux_error'
    else:
        return None

    try:
        times = np.asarray(gaia_row[time_col], dtype=float)
        fluxes = np.asarray(gaia_row[flux_col], dtype=float)
        flux_errs = np.asarray(gaia_row[flux_err_col], dtype=float)
        mags = np.asarray(gaia_row[mag_col], dtype=float)
    except (KeyError, ValueError):
        return None

    # Basic mask for valid values
    mask = ~(np.isnan(times) | np.isnan(fluxes) | np.isnan(flux_errs) | np.isnan(mags))
    if mask.sum() < min_points:
        return None

    times = times[mask]
    fluxes = fluxes[mask]
    flux_errs = flux_errs[mask]
    mags = mags[mask]

    # Convert to rest-frame time (BJD -> rest-frame days)
    times_rest = times / (1.0 + z)

    # Magnitude errors from flux errors
    mag_errs = np.where(fluxes > 0.0, K_MAG * flux_errs / fluxes, np.nan)
    bad = np.isnan(mag_errs)

    if bad.all():
        return None

    # Drop epochs with undefined mag errors
    if bad.any():
        keep = ~bad
        times_rest = times_rest[keep]
        mags = mags[keep]
        mag_errs = mag_errs[keep]

    if len(times_rest) < min_points:
        return None

    # Sort by time (celerite2 requires sorted times)
    sort_idx = np.argsort(times_rest)
    times_rest = times_rest[sort_idx]
    mags = mags[sort_idx]
    mag_errs = mag_errs[sort_idx]

    return {
        'times_rest': times_rest,
        'mags': mags,
        'mag_errs': mag_errs
    }


def _gaia_loo_process_one(i, gaia_table, catalog_table, band, min_points, max_passes):
    """
    Process a single source for Gaia LOO outlier detection.

    Returns
    -------
    tuple or None
        (result_dict, removed_orig, n_full_fail, n_loo_fail, n_post_fail, source_id)
        or None if skipped
    """
    gaia_row = gaia_table[i]
    cat_row = catalog_table[i]

    source_id = int(gaia_row['source_id'])
    z = float(cat_row['z'])

    # Prepare lightcurve
    lc = prepare_gaia_lightcurve(gaia_row, z, band=band, min_points=min_points)
    if lc is None:
        return None

    times_rest = np.asarray(lc['times_rest'], dtype=float)
    mags = np.asarray(lc['mags'], dtype=float)
    mag_errs = np.asarray(lc['mag_errs'], dtype=float)

    if len(times_rest) < min_points:
        return None

    # Center magnitudes
    mags_centered = mags - np.nanmedian(mags)

    # Arrays that will be modified as points are removed
    times_cur = times_rest.copy()
    mags_cur = mags_centered.copy()
    errs_cur = mag_errs.copy()
    orig_indices = np.arange(len(times_rest), dtype=int)

    removed_orig = []
    n_removed = 0
    n_passes_done = 0

    first_full_fit = None
    last_post_fit = None

    n_full_fail = 0
    n_loo_fail = 0
    n_post_fail = 0

    while True:
        N_cur = len(times_cur)
        if N_cur < min_points:
            break

        # 1. Fit DRW on current lightcurve
        full_fit = fit_drw_unconstrained(times_cur, mags_cur, errs_cur)
        if full_fit is None or not full_fit.get('success', True):
            n_full_fail += 1
            last_post_fit = None
            break

        # Store the very first full fit
        if first_full_fit is None:
            first_full_fit = full_fit

        log_sigma_full_cur = full_fit['log_sigma']
        log_tau_full_cur = full_fit['log_tau']
        logD_full_cur = 2.0 * log_sigma_full_cur - log_tau_full_cur

        # 2. Leave-one-out fits
        delta_logD = np.full(N_cur, np.nan, dtype=float)
        n_failed_subfits = 0

        for j in range(N_cur):
            mask = np.ones(N_cur, dtype=bool)
            mask[j] = False

            fit_j = fit_drw_unconstrained(times_cur[mask], mags_cur[mask], errs_cur[mask])
            if fit_j is None or not fit_j.get('success', True):
                n_failed_subfits += 1
                continue

            log_sigma_j = fit_j['log_sigma']
            log_tau_j = fit_j['log_tau']
            logD_j = 2.0 * log_sigma_j - log_tau_j
            delta_logD[j] = logD_j - logD_full_cur

        frac_failed = n_failed_subfits / float(N_cur)
        if frac_failed > MAX_FAILED_FRACTION or np.all(~np.isfinite(delta_logD)):
            n_loo_fail += 1
            last_post_fit = None
            break

        idx_max_D_cur = int(np.nanargmax(np.abs(delta_logD)))
        delta_best = float(delta_logD[idx_max_D_cur])

        # Decide whether to remove this point
        remove_point = True
        if ONLY_IF_D_DECREASES and not (delta_best < 0.0):
            remove_point = False
        if np.abs(delta_best) < MIN_ABS_DELTA_LOGD:
            remove_point = False

        if not remove_point:
            last_post_fit = full_fit
            break

        # 3. Refit after removing most influential point
        mask_keep = np.ones(N_cur, dtype=bool)
        mask_keep[idx_max_D_cur] = False

        post_fit = fit_drw_unconstrained(times_cur[mask_keep],
                                         mags_cur[mask_keep],
                                         errs_cur[mask_keep])
        if post_fit is None or not post_fit.get('success', True):
            n_post_fail += 1
            last_post_fit = None
            break

        # Record removed index
        removed_orig.append(int(orig_indices[idx_max_D_cur]))

        # Update current arrays
        times_cur = times_cur[mask_keep]
        mags_cur = mags_cur[mask_keep]
        errs_cur = errs_cur[mask_keep]
        orig_indices = orig_indices[mask_keep]

        last_post_fit = post_fit
        n_removed += 1
        n_passes_done += 1

        if n_passes_done >= max_passes:
            break

    # Skip if no valid fits
    if first_full_fit is None:
        return None

    # If LOO never accepted any removal, use full fit as post fit
    if last_post_fit is None:
        last_post_fit = first_full_fit

    log_sigma_full = first_full_fit['log_sigma']
    log_tau_full = first_full_fit['log_tau']
    logD_full = 2.0 * log_sigma_full - log_tau_full

    log_sigma_post = last_post_fit['log_sigma']
    log_tau_post = last_post_fit['log_tau']
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


def run_gaia_loo_for_band(gaia_table, catalog_table, band='g',
                          min_points=MIN_POINTS, max_passes=MAX_LOO_PASSES,
                          n_jobs=-1):
    """
    Run LOO outlier detection for one Gaia band using parallel processing.

    Parameters
    ----------
    gaia_table : astropy.table.Table
        Gaia lightcurve table
    catalog_table : astropy.table.Table
        Catalog table with z column
    band : str
        Band to process ('g', 'bp', or 'rp')
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

    n_obj = len(gaia_table)
    indices = np.arange(n_obj)

    results = []
    removed_indices = {}
    total_full_fail = 0
    total_loo_fail = 0
    total_post_fail = 0

    def wrapped(i):
        return _gaia_loo_process_one(i, gaia_table, catalog_table, band, min_points, max_passes)

    print(f"Running Gaia LOO for band {band} with joblib on {n_obj} objects, n_jobs={n_jobs}")

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

    print(f"Finished Gaia LOO for band {band}")
    print(f"  Successful fits: {len(results)}/{n_obj}")
    print(f"  Full-fit failures: {total_full_fail}")
    print(f"  LOO failures: {total_loo_fail}")
    print(f"  Post-fit failures: {total_post_fail}")

    results_table = Table(rows=results)

    return results_table, removed_indices


def create_loo_cleaned_lightcurves(gaia_table, removed_dict, band):
    """
    Create cleaned Gaia lightcurve table with LOO-identified outliers removed.

    Parameters
    ----------
    gaia_table : astropy.table.Table
        Original Gaia lightcurve table
    removed_dict : dict
        Dictionary mapping source_id -> array of removed indices
    band : str
        Band to clean ('g', 'bp', or 'rp')

    Returns
    -------
    gaia_cleaned : astropy.table.Table
        Cleaned lightcurve table
    """
    print("\nCreating LOO-cleaned lightcurve table...")
    print(f"  Processing {band}-band...")

    # Map band to column names
    if band == 'g':
        time_col = 'g_transit_time'
        mag_col = 'g_transit_mag'
        flux_col = 'g_transit_flux'
        flux_err_col = 'g_transit_flux_error'
    elif band == 'bp':
        time_col = 'bp_obs_time'
        mag_col = 'bp_mag'
        flux_col = 'bp_flux'
        flux_err_col = 'bp_flux_error'
    elif band == 'rp':
        time_col = 'rp_obs_time'
        mag_col = 'rp_mag'
        flux_col = 'rp_flux'
        flux_err_col = 'rp_flux_error'

    # Create copy
    gaia_cleaned = gaia_table.copy()

    n_removed_total = 0

    for i in range(len(gaia_cleaned)):
        source_id = int(gaia_cleaned['source_id'][i])

        if source_id not in removed_dict:
            continue

        removed_indices = removed_dict[source_id]
        if len(removed_indices) == 0:
            continue

        # Create mask to keep all except removed indices
        n_epochs = len(gaia_cleaned[time_col][i])
        keep_mask = np.ones(n_epochs, dtype=bool)
        keep_mask[removed_indices] = False

        n_removed_total += len(removed_indices)

        # Apply mask to all columns
        gaia_cleaned[time_col][i] = np.array(gaia_cleaned[time_col][i][keep_mask], dtype=float)
        gaia_cleaned[mag_col][i] = np.array(gaia_cleaned[mag_col][i][keep_mask], dtype=float)
        gaia_cleaned[flux_col][i] = np.array(gaia_cleaned[flux_col][i][keep_mask], dtype=float)
        gaia_cleaned[flux_err_col][i] = np.array(gaia_cleaned[flux_err_col][i][keep_mask], dtype=float)

    print(f"    {band}-band: removed {n_removed_total} epochs total")

    return gaia_cleaned


def main():
    parser = argparse.ArgumentParser(
        description="Apply LOO outlier cleaning to Gaia lightcurves",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-gaia', type=str, default='data/gaia_ztf_qso_sample_gaia_lc.fits',
                       help='Input Gaia lightcurves')
    parser.add_argument('--input-catalog', type=str, default='data/gaia_ztf_qso_sample_catalog.fits',
                       help='Input catalog (for z)')
    parser.add_argument('--output-results', type=str, default=None,
                       help='Output LOO results table (auto-generated if not provided)')
    parser.add_argument('--output-cleaned', type=str, default=None,
                       help='Output LOO-cleaned lightcurves (auto-generated if not provided)')
    parser.add_argument('--band', type=str, required=True, choices=['g', 'bp', 'rp'],
                       help='Band to process (g, bp, or rp) - REQUIRED')
    parser.add_argument('--n-jobs', type=int, default=-1, help='Number of parallel jobs (-1=all)')
    parser.add_argument('--test', action='store_true', help='Test mode (use subset)')
    parser.add_argument('--n-test', type=int, default=100, help='Number of sources for test mode')

    args = parser.parse_args()

    # Single band processing
    band = args.band.lower()

    # Auto-generate output filenames if not provided
    if args.output_results is None:
        args.output_results = f'data/gaia_loo_results_{band}.fits'
    if args.output_cleaned is None:
        args.output_cleaned = f'data/gaia_ztf_qso_sample_gaia_lc_{band}_clean_loo.fits'

    print("="*70)
    print("GAIA LOO OUTLIER CLEANING (SINGLE BAND)")
    print("="*70)
    print(f"Input Gaia:    {args.input_gaia}")
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
    gaia = Table.read(args.input_gaia)
    catalog = Table.read(args.input_catalog)
    print(f"  Gaia: {len(gaia)} sources")
    print(f"  Catalog: {len(catalog)} sources")

    # Verify alignment
    assert np.all(gaia['source_id'] == catalog['source_id']), "Gaia-Catalog alignment failed!"
    print("  ✓ Gaia-Catalog aligned")

    # Test mode subset
    if args.test:
        gaia = gaia[:args.n_test]
        catalog = catalog[:args.n_test]
        print(f"\nTest mode: processing {len(gaia)} sources")

    # Run LOO for single band
    print(f"\n{'='*70}")
    print(f"Processing {band}-band")
    print("="*70)

    results_table, removed_dict = run_gaia_loo_for_band(
        gaia, catalog, band=band, n_jobs=args.n_jobs
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
    gaia_cleaned = create_loo_cleaned_lightcurves(gaia, removed_dict, band)

    gaia_cleaned.write(args.output_cleaned, format='fits', overwrite=True)
    file_size_mb = os.path.getsize(args.output_cleaned) / (1024 * 1024)
    print(f"\n✓ Saved LOO-cleaned lightcurves: {args.output_cleaned}")
    print(f"  {len(gaia_cleaned)} sources ({file_size_mb:.1f} MB)")

    print("\n" + "="*70)
    print("LOO CLEANING COMPLETE")
    print("="*70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
