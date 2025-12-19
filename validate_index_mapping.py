#!/usr/bin/env python
"""
Validate LOO index mapping fix.

Tests that sorted_to_orig mapping correctly identifies outliers in
original FITS array positions after sorting.

Usage:
    # Unit test
    python validate_index_mapping.py --unit-test

    # Validate actual cleaned lightcurves
    python validate_index_mapping.py \
        --ztf-original data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits \
        --ztf-cleaned data/gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits \
        --loo-results data/ztf_loo_results_g.fits \
        --band g \
        --n-test 100
"""

import sys
import numpy as np
import argparse
from astropy.table import Table

from apply_loo_cleaning import prepare_ztf_lightcurve


def test_unit_index_mapping():
    """
    Unit test: Verify sorted_to_orig mapping works correctly.
    """
    print("="*70)
    print("UNIT TEST: Index Mapping")
    print("="*70)

    # Create mock FITS row with UNSORTED times
    class MockRow:
        def __init__(self, times, mags, magerrs):
            # Store as numpy arrays (matching astropy.table behavior)
            self.mjd_g = np.array(times, dtype=float)
            self.mag_g = np.array(mags, dtype=float)
            self.magerr_g = np.array(magerrs, dtype=float)

        def __getitem__(self, key):
            return getattr(self, key)

    # Test case 1: Simple unsorted array
    print("\nTest 1: Simple unsorted array")
    times = [3.0, 1.0, 5.0, 2.0, 4.0]
    mags = [19.1, 19.3, 19.0, 19.2, 19.1]
    magerrs = [0.1, 0.1, 0.1, 0.1, 0.1]

    row = MockRow(times, mags, magerrs)
    lc = prepare_ztf_lightcurve(row, z=0.5, band='g', min_points=5)

    # Verify sorted times
    times_sorted = lc['times_rest'] * 1.5  # Undo rest-frame conversion
    expected_sorted = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert np.allclose(times_sorted, expected_sorted), "Times not sorted correctly!"
    print(f"  ✓ Times sorted: {times_sorted}")

    # Verify mapping
    # Original: [3.0, 1.0, 5.0, 2.0, 4.0] (indices: 0, 1, 2, 3, 4)
    # Sorted:   [1.0, 2.0, 3.0, 4.0, 5.0] (from orig: 1, 3, 0, 4, 2)
    # So sorted_to_orig = [1, 3, 0, 4, 2]
    expected_mapping = np.array([1, 3, 0, 4, 2])
    assert np.array_equal(lc['sorted_to_orig'], expected_mapping), \
        f"Mapping wrong! Got {lc['sorted_to_orig']}, expected {expected_mapping}"
    print(f"  ✓ Mapping correct: {lc['sorted_to_orig']}")

    # Test case 2: Array with outlier at beginning
    print("\nTest 2: Outlier at beginning of original array")
    times = [100.0, 1.0, 2.0, 3.0, 4.0]  # Outlier at index 0
    mags = [18.0, 19.1, 19.2, 19.1, 19.0]
    magerrs = [0.1, 0.1, 0.1, 0.1, 0.1]

    row = MockRow(times, mags, magerrs)
    lc = prepare_ztf_lightcurve(row, z=0.0, band='g', min_points=5)

    # Sorted: [1.0, 2.0, 3.0, 4.0, 100.0]
    # From orig: [1, 2, 3, 4, 0]
    expected_mapping = np.array([1, 2, 3, 4, 0])
    assert np.array_equal(lc['sorted_to_orig'], expected_mapping), \
        f"Mapping wrong! Got {lc['sorted_to_orig']}, expected {expected_mapping}"
    print(f"  ✓ Mapping correct: {lc['sorted_to_orig']}")
    print(f"  ✓ Outlier at sorted index 4 maps to original index {lc['sorted_to_orig'][4]}")

    # Test case 3: Array with NaN values
    print("\nTest 3: Array with NaN values")
    times = [3.0, np.nan, 1.0, 5.0, 2.0]
    mags = [19.1, 19.3, 19.3, 19.0, 19.2]
    magerrs = [0.1, 0.1, 0.1, 0.1, 0.1]

    row = MockRow(times, mags, magerrs)
    lc = prepare_ztf_lightcurve(row, z=0.0, band='g', min_points=4)

    # After NaN filter: indices [0, 2, 3, 4] survive
    # Times: [3.0, 1.0, 5.0, 2.0]
    # Sorted: [1.0, 2.0, 3.0, 5.0]
    # From filtered: [2, 4, 0, 3] (indices in filtered array)
    # Mapped to original: [2, 4, 0, 3]
    expected_mapping = np.array([2, 4, 0, 3])
    assert np.array_equal(lc['sorted_to_orig'], expected_mapping), \
        f"Mapping wrong! Got {lc['sorted_to_orig']}, expected {expected_mapping}"
    print(f"  ✓ Mapping correct (with NaN): {lc['sorted_to_orig']}")

    print("\n" + "="*70)
    print("✅ ALL UNIT TESTS PASSED")
    print("="*70)


def validate_cleaned_lightcurves(ztf_original, ztf_cleaned, loo_results, band='g', n_test=100):
    """
    Validate that cleaned lightcurves have correct epochs removed.

    Parameters
    ----------
    ztf_original : astropy.table.Table
        Original ZTF lightcurves (before LOO cleaning)
    ztf_cleaned : astropy.table.Table
        Cleaned ZTF lightcurves (after LOO cleaning)
    loo_results : astropy.table.Table
        LOO results with n_removed column
    band : str
        Band to validate ('g' or 'r')
    n_test : int
        Number of sources to test
    """
    print("\n" + "="*70)
    print(f"VALIDATION: Cleaned Lightcurves ({band}-band)")
    print("="*70)

    n_test = min(n_test, len(ztf_original), len(ztf_cleaned), len(loo_results))
    print(f"Testing {n_test} sources...")

    n_pass = 0
    n_fail = 0
    issues = []

    for i in range(n_test):
        source_id_orig = ztf_original['source_id'][i]
        source_id_clean = ztf_cleaned['source_id'][i]

        # Verify alignment
        if source_id_orig != source_id_clean:
            issues.append(f"Source {i}: ID mismatch ({source_id_orig} vs {source_id_clean})")
            n_fail += 1
            continue

        # Get LOO info
        loo_row = loo_results[loo_results['source_id'] == source_id_orig]
        if len(loo_row) == 0:
            continue  # No LOO result for this source

        n_removed = int(loo_row['n_removed'][0])

        # Check epoch counts
        n_orig = len(ztf_original[f'mjd_{band}'][i])
        n_clean = len(ztf_cleaned[f'mjd_{band}'][i])
        expected_clean = n_orig - n_removed

        if n_clean != expected_clean:
            issues.append(f"Source {source_id_orig}: Expected {expected_clean} epochs, got {n_clean}")
            n_fail += 1
        else:
            n_pass += 1

        # Progress
        if (i + 1) % 20 == 0:
            print(f"  Tested {i+1}/{n_test} sources...")

    print(f"\nResults:")
    print(f"  ✅ Passed: {n_pass}/{n_test}")
    print(f"  ❌ Failed: {n_fail}/{n_test}")

    if n_fail > 0:
        print(f"\nFirst 10 issues:")
        for issue in issues[:10]:
            print(f"  - {issue}")

    if n_fail == 0:
        print("\n" + "="*70)
        print("✅ VALIDATION PASSED: All cleaned lightcurves have correct epoch counts")
        print("="*70)
    else:
        print("\n" + "="*70)
        print("❌ VALIDATION FAILED: Some cleaned lightcurves have incorrect epoch counts")
        print("="*70)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Validate LOO index mapping fix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--unit-test', action='store_true',
                       help='Run unit tests only')
    parser.add_argument('--ztf-original', type=str,
                       help='Original ZTF lightcurves FITS file')
    parser.add_argument('--ztf-cleaned', type=str,
                       help='Cleaned ZTF lightcurves FITS file')
    parser.add_argument('--loo-results', type=str,
                       help='LOO results FITS file')
    parser.add_argument('--band', type=str, default='g', choices=['g', 'r'],
                       help='Band to validate')
    parser.add_argument('--n-test', type=int, default=100,
                       help='Number of sources to test')

    args = parser.parse_args()

    if args.unit_test:
        test_unit_index_mapping()
        return 0

    if not (args.ztf_original and args.ztf_cleaned and args.loo_results):
        print("Error: --ztf-original, --ztf-cleaned, and --loo-results required")
        print("       (or use --unit-test for unit tests only)")
        sys.exit(1)

    # Run unit tests first
    test_unit_index_mapping()

    # Load data
    print("\nLoading data...")
    ztf_orig = Table.read(args.ztf_original)
    ztf_clean = Table.read(args.ztf_cleaned)
    loo_res = Table.read(args.loo_results)
    print(f"  Original: {len(ztf_orig)} sources")
    print(f"  Cleaned:  {len(ztf_clean)} sources")
    print(f"  LOO:      {len(loo_res)} sources")

    # Validate
    validate_cleaned_lightcurves(ztf_orig, ztf_clean, loo_res, args.band, args.n_test)

    return 0


if __name__ == '__main__':
    sys.exit(main())
