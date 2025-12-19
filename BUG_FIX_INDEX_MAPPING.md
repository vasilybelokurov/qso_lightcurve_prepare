# Bug Fix: LOO Outlier Index Mapping (2025-12-10)

## Executive Summary

**Critical bug fixed:** LOO outlier removal code was applying indices from **sorted** arrays to **unsorted** original FITS arrays, resulting in removal of **incorrect epochs**.

**Impact:**
- ✅ LOO analysis and DRW fits: **CORRECT** (operated on sorted data)
- ❌ Cleaned lightcurve exports: **WRONG EPOCHS REMOVED**
- ❌ Subsequent re-fits with different parameters: **CORRUPTED**

**Status:** Fixed in `apply_loo_cleaning.py` (2025-12-10)

---

## Root Cause

### The Problem

**ZTF lightcurves are NOT pre-sorted by time** in FITS files. Arrays contain time jumps:

```python
# Example from data/matched_5k_ztf_lc.fits, source 5:
MJD: [58298, 58301, ..., 58347, 60610, 58347, ...]
                                     ↑ Jump from 58347 → 60610 → 58347
```

**The code workflow:**

1. `prepare_ztf_lightcurve()` **sorts** arrays (line 84-88):
   ```python
   sort_idx = np.argsort(times)
   times = times[sort_idx]  # SORTED
   ```

2. LOO analysis operates on **sorted** arrays (lines 132-234):
   ```python
   orig_indices = np.arange(len(times_rest))  # Indices into SORTED array
   removed_orig.append(orig_indices[idx_max_D_cur])  # BUG: sorted index!
   ```

3. `create_loo_cleaned_lightcurves()` applies to **unsorted** FITS arrays (lines 385-399):
   ```python
   removed_idx = removed_dict[source_id]  # Sorted indices!
   keep_mask[removed_idx] = False  # Applied to UNSORTED array!
   mjd_clean.append(ztf_table['mjd_g'][i][keep_mask])  # WRONG EPOCHS!
   ```

---

## Visual Example

**Original UNSORTED FITS array:**
```
Index:  0     1     2     3     4     5
MJD:    58300 58290 58310 58305 58295 60610
Mag:    19.1  19.3  19.0  19.8  19.2  18.5  ← Outlier (60610 = index 5)
```

**After prepare_ztf_lightcurve sorting:**
```
Sorted:       0     1     2     3     4     5
MJD:          58290 58295 58300 58305 58310 60610
Mag:          19.3  19.2  19.1  19.8  19.0  18.5
Orig Index:   1     4     0     3     2     5     ← Mapping
```

**LOO identifies sorted index 5 as outlier** ✅ Correct (MJD=60610)

**BUGGY behavior (before fix):**
```python
removed_indices[source_id] = [5]  # Index in SORTED array

# Applied to UNSORTED array:
keep_mask[5] = False  # Removes position 5 in UNSORTED
# Position 5 in UNSORTED = MJD 60610 ✅ Happens to be correct (by luck!)

# BUT if original order was different:
# UNSORTED: [60610, 58290, 58295, 58300, 58305, 58310]
# keep_mask[5] = False removes MJD 58310 ❌ WRONG!
```

**The bug manifests when original array order ≠ sorted order at the outlier position.**

---

## The Fix

### Changed Files

#### 1. `apply_loo_cleaning.py::prepare_ztf_lightcurve()` (Lines 40-112)

**Added:**
- Track original FITS indices through NaN filtering: `orig_indices_after_mask`
- Map sorted indices back to original: `sorted_to_orig`
- Return mapping in dictionary: `'sorted_to_orig': sorted_to_orig`

**Key code:**
```python
# Line 88: Track which indices survived NaN filter
orig_indices_after_mask = np.where(mask)[0]

# Line 95-102: Create mapping from sorted to original
sort_idx = np.argsort(times)
times = times[sort_idx]
sorted_to_orig = orig_indices_after_mask[sort_idx]

# Line 111: Return mapping
return {'times_rest': times_rest, ..., 'sorted_to_orig': sorted_to_orig}
```

#### 2. `apply_loo_cleaning.py::_ztf_loo_process_one()` (Lines 115-269)

**Changed:**
- Receive `sorted_to_orig` from `prepare_ztf_lightcurve()` (line 138)
- Track positions in sorted array: `sorted_indices` (line 149)
- Map removed indices to original FITS positions (lines 226-229)

**Key code:**
```python
# Line 138: Get mapping
sorted_to_orig = lc["sorted_to_orig"]

# Line 226-229: Map sorted index → original FITS index
sorted_idx_to_remove = sorted_indices[idx_max_D_cur]
orig_idx_to_remove = sorted_to_orig[sorted_idx_to_remove]
removed_orig.append(int(orig_idx_to_remove))
```

#### 3. `apply_loo_cleaning.py::create_loo_cleaned_lightcurves()` (Lines 351-433)

**No code changes** - now receives **correct** original FITS indices!

**Added documentation:**
- Updated docstring to clarify indices refer to original FITS arrays (lines 360-364)
- Added note about bug fix date (lines 373-377)

#### 4. `fit_gaia_loo_parallel_bfgs.py` (Lines 1-170)

**Added defensive programming:**
- Docstring warning about sorting assumption (lines 14-20)
- Runtime check that Gaia arrays are sorted (lines 158-163)

**Note:** Gaia DR3 FITS arrays **are pre-sorted**, so this code is safe. Check added for future-proofing.

---

## Validation

### Test Case: Synthetic Example

**Setup:**
```python
# Create test lightcurve with known outlier
times_unsorted = np.array([3.0, 1.0, 5.0, 2.0, 4.0, 100.0])  # Outlier at index 5
mags = np.array([19.1, 19.3, 19.0, 19.2, 19.1, 18.0])

# After sorting:
# times_sorted = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0]
# sorted_to_orig = [1, 3, 0, 4, 2, 5]
# LOO identifies sorted index 5 (time=100.0) as outlier
```

**Expected:**
```python
removed_indices = [5]  # Original FITS index
times_clean = np.delete(times_unsorted, 5)  # Removes time=100.0 ✅
```

**Verification script:** `validate_index_mapping.py` (see below)

---

## Affected Output Files

### Corrupted Files (Created Before Fix)

**If you ran `apply_loo_cleaning.py` before 2025-12-10:**

| File | Status | Action Required |
|------|--------|-----------------|
| `data/gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits` | ❌ **CORRUPTED** | **Regenerate** |
| `data/gaia_ztf_qso_sample_ztf_lc_r_clean_loo.fits` | ❌ **CORRUPTED** | **Regenerate** |
| `data/ztf_loo_results_g.fits` | ✅ **VALID** | DRW parameters correct |
| `data/ztf_loo_results_r.fits` | ✅ **VALID** | DRW parameters correct |

### Not Affected

**These files are SAFE (Gaia data pre-sorted):**
- `data/gaia_loo_full_results_*.npz` ✅
- All Gaia LOO results ✅
- All analysis plots ✅ (used parameters, not cleaned lightcurves)

---

## Regeneration Instructions

### Step 1: Re-run LOO with Fixed Code

```bash
# Activate environment
source ~/Work/venvs/.venv/bin/activate

# Process g-band
cd ~/IoA\ Dropbox/Dr\ V.A.\ Belokurov/Code/qso_lightcurve_prepare
python apply_loo_cleaning.py \
    --input-ztf data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits \
    --input-catalog data/gaia_ztf_qso_sample_catalog.fits \
    --band g \
    --n-jobs 8

# Process r-band
python apply_loo_cleaning.py \
    --input-ztf data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits \
    --input-catalog data/gaia_ztf_qso_sample_catalog.fits \
    --band r \
    --n-jobs 8
```

### Step 2: Verify Cleaned Output

```bash
# Run validation script
python validate_index_mapping.py \
    --cleaned-ztf data/gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits \
    --removed-indices data/ztf_loo_results_g.fits
```

### Step 3: Re-run Downstream Analysis

**Any analysis using cleaned lightcurves must be re-run:**
- Re-fits with different jitter values
- Parameter recovery tests
- Comparison plots using cleaned data

---

## Technical Details

### Index Mapping Chain

```
Original FITS indices (N epochs)
    ↓ NaN filter (mask)
Filtered indices (M epochs, M ≤ N)
    ↓ Time sort (sort_idx)
Sorted indices (M epochs, reordered)
    ↓ LOO analysis
Identified outliers (sorted positions)
    ↓ sorted_to_orig mapping [THE FIX]
Original FITS indices of outliers
    ↓ Apply to FITS arrays
Cleaned lightcurves ✅
```

### Example: Full Mapping

**Original FITS (N=6):**
```
Index:  0     1     2     3     4     5
MJD:    58300 58290 58310 58305 58295 60610
Valid:  ✓     ✓     ✓     ✓     ✓     ✓     (no NaNs)
```

**After NaN filter (M=6, no change):**
```
orig_indices_after_mask = [0, 1, 2, 3, 4, 5]
```

**After sorting:**
```
Sorted:         0     1     2     3     4     5
MJD:            58290 58295 58300 58305 58310 60610
sort_idx:       1     4     0     3     2     5
sorted_to_orig: 1     4     0     3     2     5     ← THE KEY MAPPING
```

**LOO identifies sorted index 5 as outlier:**
```
sorted_idx_to_remove = 5
orig_idx_to_remove = sorted_to_orig[5] = 5  ✅ Correct!
```

**Edge case (different original order):**
```
Original:  [60610, 58290, 58295, 58300, 58305, 58310]
Sorted:    [58290, 58295, 58300, 58305, 58310, 60610]
sort_idx:  [1,     2,     3,     4,     5,     0    ]
sorted_to_orig: [1, 2, 3, 4, 5, 0]

LOO identifies sorted index 5 (MJD=60610):
sorted_idx_to_remove = 5
orig_idx_to_remove = sorted_to_orig[5] = 0  ✅ Correct! (60610 at position 0)
```

---

## Testing

### Unit Test

Create `test_index_mapping.py`:

```python
import numpy as np
from apply_loo_cleaning import prepare_ztf_lightcurve

def test_index_mapping():
    """Test that sorted_to_orig correctly maps indices."""
    # Create mock FITS row with unsorted times
    class MockRow:
        def __init__(self):
            self.data = {
                'mjd_g': np.array([3.0, 1.0, 5.0, 2.0, 4.0]),
                'mag_g': np.array([19.1, 19.3, 19.0, 19.2, 19.1]),
                'magerr_g': np.array([0.1, 0.1, 0.1, 0.1, 0.1]),
            }
        def __getitem__(self, key):
            return self.data[key]

    row = MockRow()
    lc = prepare_ztf_lightcurve(row, z=0.5, band='g')

    # Verify sorted times
    expected_times_sorted = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    np.testing.assert_array_equal(lc['times_rest'] * 1.5, expected_times_sorted)

    # Verify mapping: sorted_to_orig[sorted_idx] = original_idx
    # Sorted:  [1.0, 2.0, 3.0, 4.0, 5.0]
    # Original: [3.0, 1.0, 5.0, 2.0, 4.0]
    # Mapping:  [1,   3,   0,   4,   2  ]
    expected_mapping = np.array([1, 3, 0, 4, 2])
    np.testing.assert_array_equal(lc['sorted_to_orig'], expected_mapping)

    print("✅ Index mapping test PASSED")

if __name__ == '__main__':
    test_index_mapping()
```

**Run:** `python test_index_mapping.py`

---

## References

- Original bug report: User message, 2025-12-10
- ZTF data verification: `data/matched_5k_ztf_lc.fits` source 5 (MJD jump at index 453)
- Gaia pre-sorted verification: All 10 tested sources in `gaia/quaia_gaia_lightcurves.fits`

---

## Changelog

**2025-12-10:** Initial bug fix
- Fixed `prepare_ztf_lightcurve()` to return `sorted_to_orig` mapping
- Fixed `_ztf_loo_process_one()` to map sorted indices to original FITS indices
- Added defensive check to `fit_gaia_loo_parallel_bfgs.py`
- Created this documentation

---

**Status:** ✅ **FIXED** - Code correct as of 2025-12-10
**Action Required:** Regenerate all cleaned ZTF lightcurve FITS files
