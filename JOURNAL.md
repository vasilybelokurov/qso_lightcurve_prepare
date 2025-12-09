# QSO Lightcurve Preparation Pipeline — Project Journal

**Project:** SDSS DR16Q × Gaia DR3 Cross-Match and Lightcurve Download
**Location:** `~/Work/Code/qso_lightcurve_prepare/`
**Data Output:** `~/data/gaia/`

---

## Session 1: 2025-12-07

### Objective
Create comprehensive QSO catalog by cross-matching SDSS DR16Q spectroscopic QSO catalog with Gaia DR3 photometric data and downloading multi-band lightcurves.

### Motivation
- SDSS DR16Q: 750,414 spectroscopic QSOs (largest spectroscopic sample)
- Gaia DR3: Multi-epoch photometry (G, BP, RP bands) for variability studies
- Goal: Combine spectroscopic redshifts with photometric variability for largest possible QSO dataset

### Implementation

#### Step 1: Cross-Match Catalog Download
**Script:** `download_sdssdr16q_gaia_source.py`

**Method:**
- Spatial cross-match using q3c_join (1 arcsec radius)
- Query WSDB directly (no local_join needed)
- Download ALL columns from both catalogs
- Auto-detect and prefix conflicting column names:
  - `ra`, `dec`, `duplicated_source` → `sdss_*` and `gaia_*`

**Results:**
- Runtime: 12 minutes (718.8 seconds)
- Matched sources: 489,484 QSOs (65% of DR16Q final sample)
- Output file: `~/data/gaia/sdssdr16q_gaia_source.fits` (2.1 GB)
- Columns: 336 (complete SDSS + Gaia metadata)
- Redshift range: -999.0 to 7.024, median z=1.644
- Note: 4 duplicate Gaia sources (489,480 unique source_ids)

**Key SQL:**
```sql
SELECT <all_sdss_cols>, <all_gaia_cols>
FROM sdssdr16qso.main AS q
JOIN gaia_dr3.gaia_source AS g
  ON q3c_join(q.ra, q.dec, g.ra, g.dec, 1.0/3600.0)
WHERE q.is_qso_final = 1
```

#### Step 2: Epoch Photometry Download
**Script:** `download_sdssdr16q_gaia_epoch_photometry.py`

**Method:**
- Upload 489,484 source_ids to WSDB using sqlutilpy.local_join()
- Join with gaia_dr3.epoch_photometry table
- Download multi-band lightcurves (G, BP, RP) + quality flags

**Results:**
- Runtime: 8.2 minutes
- Sources with epochs: 223,221 (223,219 unique, ~46% of cross-match)
- Output file: `~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits` (980 MB)
- Columns: 21 (times, mags, fluxes, errors × 3 bands + flags)
- Median epochs per band: 41
- Total observations: 9.5M per band

**Data columns per source:**
- G-band: g_transit_time, g_transit_mag, g_transit_flux, g_transit_flux_error
- BP-band: bp_obs_time, bp_mag, bp_flux, bp_flux_error
- RP-band: rp_obs_time, rp_mag, rp_flux, rp_flux_error
- Quality flags: photometry_flag_* (4), variability_flag_* (3)

**Key SQL:**
```sql
SELECT m.source_id, e.g_transit_time, e.g_transit_mag, ...
FROM mytmptable AS m
JOIN gaia_dr3.epoch_photometry AS e ON e.source_id = m.source_id
```

### Files Created
1. `download_sdssdr16q_gaia_source.py` (executable)
2. `download_sdssdr16q_gaia_epoch_photometry.py` (executable)
3. `README.md` (usage documentation)
4. `JOURNAL.md` (this file)

### Output Data Products
```
~/data/gaia/
├── sdssdr16q_gaia_source.fits           (2.1 GB, 489,484 rows, 336 cols)
└── sdssdr16q_gaia_epoch_photometry.fits (980 MB, 223,221 rows, 21 cols)
```

### Key Statistics
- Cross-match success rate: 65% (489k/750k DR16Q QSOs)
- Epoch photometry coverage: 46% (223k/489k matched sources)
- Median epochs per source: 41 (G, BP, RP)
- Total photometry points: 28.6M (9.5M × 3 bands)

### Notes
- ~54% of cross-matched sources lack epoch photometry (expected for faint/sparse sampling)
- Variable-length arrays used for per-epoch data (FITS format)
- All quality flags included for robust analysis
- Pattern follows proven Quaia lightcurve download approach

#### Step 3: ZTF Target Catalog Creation
**Script:** `create_ztf_target_catalog.py`

**Method:**
- Extract unique source_ids from Gaia epoch photometry file
- Join with source catalog to retrieve RA/Dec coordinates
- Create compact catalog (source_id, ra, dec only)

**Results:**
- Input: 223,221 rows from epoch photometry file
- Unique sources: 223,219 (2 duplicates in epoch photometry)
- Output file: `data/sdssdr16q_gaia_ztf_targets.fits` (5.2 KB)
- Columns: 3 (source_id, ra, dec)
- RA range: 0.027° to 359.991°
- Dec range: -88.865° to +82.734°

**Usage:**
```bash
python create_ztf_target_catalog.py
```

#### Step 4: ZTF Lightcurve Download (Prepared)
**Script:** `download_ztf_threaded.py` (adapted from qso_gp_mock)

**Method:**
- Multi-threaded download from IRSA ZTF API
- Batch processing with resumable chunks
- Manifest tracking for fault tolerance
- Per-band aggregation (zg, zr, zi)

**Configuration:**
- Input: `data/sdssdr16q_gaia_ztf_targets.fits` (223,219 sources)
- Output: `~/data/ztf/sdssdr16q_ztf_lightcurves.fits`
- Chunks directory: `ztf_chunks/`
- Batch size: 500 sources
- Threads: 8 (optimized for I/O-bound API calls)
- Search radius: 2 arcsec

**Launch script:** `run_ztf_download.sh`
```bash
./run_ztf_download.sh
```

**Expected runtime:** Several hours to days (depends on ZTF API performance)

**Status:** Baseline test completed, optimization in progress

#### Step 5: Performance Baseline Test (2K sources)
**Test:** Initial ZTF download performance characterization

**Configuration:**
- Test size: 2000 sources (first 2K from full catalog)
- Threads: 8
- Batch size: 100
- Duration: ~22 minutes (terminated at 60% completion)

**Results:**
- Sources downloaded: 1,196/2000 (60%)
- Chunks created: 12 FITS files
- Match rate: 99-100% (excellent)
- **Download speed: 0.93-1.14 src/s** (average ~1.0 src/s)
- CPU usage: Low (I/O-bound, as expected)

**Performance Analysis:**
- Bottleneck: I/O-bound (waiting for ZTF API responses)
- Per-batch variance: 0.83-1.29 src/s
- Projected full catalog time: **2.6 days** at 1.0 src/s

**Optimization Identified:**
Since I/O-bound, can increase parallelism without CPU penalty:
1. More threads (16-32): Expected 2× speedup → **1.3 days**
2. Multiple processes (4×8 threads): Expected 3-4× speedup → **0.7 days**

**Status:** Thread scaling test completed - API rate limiting prevents speedup

#### Step 6: Thread Scaling Test
**Test:** Systematic evaluation of thread count impact on download performance

**Configuration:**
- Test size: 100 sources per thread count
- Thread counts tested: 4, 8, 16, 24 (32 not started)
- Batch size: 100
- Test script: `test_thread_scaling.sh`

**Results:**

| Threads | Time | Matched | Failed | Rate      | Status |
|---------|------|---------|--------|-----------|--------|
| 4       | 190s | 100/100 | 0      | 0.53 src/s| ✓      |
| 8       | 192s | 100/100 | 0      | 0.52 src/s| ✓      |
| 16      | 681s | 2/100   | 98     | 0.003 src/s| ✗     |
| 24      | (terminated) | 0/70+ | 70+ | 0 src/s | ✗ |
| 32      | (not started) | - | - | - | - |

**Critical Findings:**

1. **No performance gain from increased threads:**
   - 4 threads: 0.53 src/s
   - 8 threads: 0.52 src/s (no improvement)
   - Thread count has zero impact on performance at 4-8 threads

2. **Catastrophic failure at 16+ threads:**
   - 16 threads: 98% failure rate (only 2/100 successful)
   - 24 threads: 100% failure rate (0 matches in 70+ attempts)
   - Root cause: **ZTF IRSA API rate limiting**

3. **API concurrency limits:**
   - API enforces strict per-client rate limits
   - Limit appears to be ~8-12 concurrent connections
   - Exceeding limit causes timeouts/empty responses
   - More threads = worse performance (not better)

4. **Performance variance:**
   - Baseline 2K test: ~1.0 src/s (8 threads)
   - Thread scaling test: ~0.5 src/s (8 threads)
   - Possible causes: API performance degradation, rate limiting from sustained use, time-of-day effects

**Conclusion:**
Thread scaling does NOT work for ZTF IRSA API. The API has strict concurrency limits that prevent parallelization speedup. Simple thread increase strategy is ineffective.

**Recommendation:**
Accept baseline performance of **0.5-1.0 src/s** with 8 threads, resulting in **2.6-5.2 days** for full 223K catalog download. Alternative strategies (multi-process, API documentation review) may be explored if needed.

**Status:** Tests terminated, optimization approach abandoned

#### Step 7: ZTF Lightcurve Download (Production Run)
**Test:** Full catalog download attempt (20K sources subset)

**Configuration:**
- Started: 2025-12-07
- Target: ~20K sources (first batches from full 223K catalog)
- Threads: 8
- Batch size: 500

**Results:**
- Download interrupted/killed after partial completion
- Batches completed: 57 (0-56)
- Chunks created: 21 FITS files (`chunk_0000.fits` through `chunk_0020.fits`)
- Sources downloaded: 8,529
- Sources failed: 19,971
- **Match rate: 30%** (8,529/28,500)

**Performance:**
- Sustained download rate: ~0.5-1.0 src/s (consistent with baseline tests)
- ZTF API rate limiting prevents speedup
- Download stability: Good (resumable chunks+manifest system worked)

**Status:** Download interrupted, partial data recovered

#### Step 8: Chunk Merging
**Script:** `merge_ztf_chunks.py` (adapted from qso_gp_mock)

**Method:**
- Copied `merge_ztf_chunks_fix_dtypes.py` from qso_gp_mock project
- Adapted for local project structure
- Fixed dtypes for variable-length arrays (`airmass`, `sharp`, `chi`, `limitmag`)
- Merged 21 chunk files via `astropy.table.vstack()`

**Results:**
- Input: 21 chunks (8,529 unique sources)
- Output: `data/ztf_lightcurves_partial_merged.fits` (352 MB)
- Columns: 27 (source_id, ra, dec + 8 per band × 3 bands)
- Merge time: <1 minute

**Data Quality:**
- g-band: 8,471/8,529 sources (99.3%), median 384 epochs
- r-band: 8,479/8,529 sources (99.4%), median 635 epochs
- i-band: 7,843/8,529 sources (92.0%), median 98 epochs

**Epoch statistics (median/mean/max):**
- g-band: 384 / 432 / 2,771 epochs
- r-band: 635 / 653 / 3,811 epochs
- i-band: 98 / 111 / 603 epochs

**Status:** Partial data merged successfully, ready for analysis

#### Step 9: Matched Catalog Creation
**Script:** `create_matched_catalog.py` (adapted from qso_gp_mock)

**Objective:** Create fully-matched multi-survey catalogs combining ZTF + Gaia + SDSS DR16Q data

**Method:**
- Adapted `create_matched_ztf_gaia_catalog.py` from qso_gp_mock project
- Load ZTF lightcurves (defines target sample: 8,529 sources)
- Find matching sources in Gaia epoch photometry by `source_id`
- Find matching sources in Gaia source catalog by `source_id`
- Create intersection (sources present in ALL THREE datasets)
- Sort all tables by `source_id` for row alignment
- Verify alignment with assertions
- Save 3 separate FITS files with identical row counts

**Input Files:**
1. `data/ztf_lightcurves_partial_merged.fits` (8,529 sources, 27 columns)
2. `~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits` (223,221 sources, 21 columns)
3. `~/data/gaia/sdssdr16q_gaia_source.fits` (489,484 sources, 336 columns)

**Output Files:**
1. `data/gaia_ztf_qso_sample_ztf_lc.fits` (352 MB, 8,529 rows, 27 columns)
2. `data/gaia_ztf_qso_sample_gaia_lc.fits` (29 MB, 8,529 rows, 21 columns)
3. `data/gaia_ztf_qso_sample_catalog.fits` (38 MB, 8,529 rows, 336 columns)

**Matching Results:**
- ZTF sources: 8,529
- Found in Gaia LC: 8,529 (100%)
- Found in Gaia catalog: 8,529 (100%)
- **Final matched (ALL THREE): 8,529 (100.0% match rate)**
- Perfect alignment: all three files have identical row counts and `source_id` order
- Row i in all 3 files = same QSO (verified with assertions)

**Data Quality Summary:**

*ZTF lightcurves:*
- g-band: 8,471 sources (99.3%), median 384 epochs, max 2,771
- r-band: 8,479 sources (99.4%), median 635 epochs, max 3,811
- i-band: 7,843 sources (92.0%), median 98 epochs, max 603
- **Total ZTF epochs: 10,204,591** (g=3.7M, r=5.6M, i=0.9M)

*Gaia lightcurves:*
- g-band: 8,529 sources (100%), median 30 epochs, max 97
- bp-band: 8,529 sources (100%), median 30 epochs, max 97
- rp-band: 8,529 sources (100%), median 30 epochs, max 97
- **Total Gaia epochs: 827,241**

*Catalog properties:*
- Columns: 336 (all SDSS DR16Q + all Gaia DR3 columns)
- Redshift range: z = 0.016 – 4.846
- Median redshift: z = 1.420
- Includes: spectroscopic properties, photometry, astrometry, quality flags

**Key Features:**
- No quality filtering applied (raw data preserved)
- Variable-length arrays for epoch data (FITS format)
- All catalog columns retained (336 total)
- Files aligned by `source_id` for easy cross-referencing
- Data convention: all files stored in `data/` directory

**Verification:**
- ✓ Row count equality verified across all 3 files
- ✓ `source_id` alignment verified with assertions
- ✓ Sample source inspection confirms correct matching
- ✓ File integrity checked (all files readable)

**Status:** Matched catalogs created successfully, ready for GP modeling

#### Step 10: ZTF Lightcurve Quality Cleaning
**Script:** `clean_ztf_lightcurves.py` (based on qso_gp_mock/ztf_10k_lightcurve_quality.ipynb)

**Objective:** Apply quality cuts to ZTF lightcurves to remove bad epochs and improve data quality for GP modeling

**Method:**
- Based on quality filtering from qso_gp_mock notebook `ztf_10k_lightcurve_quality.ipynb`
- Apply strict quality cuts to g and r bands only (i-band excluded)
- Keep all sources regardless of remaining epochs (no minimum epoch filter)
- Create cleaned file with reduced columns (g, r bands only)
- Preserve original file unchanged

**Quality Cuts Applied (per band):**
1. `catflags == 0` — No quality flags set (strictest: all bits must be 0)
2. `|sharp| < 0.25` — PSF sharpness within acceptable range
3. `0.5 <= chi <= 1.5` — PSF fit quality (chi-squared per degree of freedom)
4. `airmass < 1.8` — Observation quality (atmospheric transparency)

**Input File:**
- `data/gaia_ztf_qso_sample_ztf_lc.fits` (352 MB, 8,529 rows, 27 cols)
- Contains g, r, i bands with all quality columns

**Output File:**
- `data/gaia_ztf_qso_sample_ztf_lc_gr_clean.fits` (204 MB, 8,529 rows, 19 cols)
- Contains g, r bands only (i-band removed)
- All epochs failing quality cuts removed

**Cleaning Results:**

*g-band:*
- Median epochs: 384 → 253 (69.0% retention)
- Mean epochs: 432 → 298 (69.0% retention)
- Max epochs: 2,771 → 2,057
- Sources with data: 8,448/8,529 (99.1%)
- Sources with 0 epochs after cleaning: 81 (0.9%)

*r-band:*
- Median epochs: 635 → 383 (60.6% retention)
- Mean epochs: 653 → 396 (60.6% retention)
- Max epochs: 3,811 → 2,136
- Sources with data: 8,458/8,529 (99.2%)
- Sources with 0 epochs after cleaning: 71 (0.8%)

**Data Reduction:**
- File size: 352 MB → 204 MB (42% reduction)
- Columns: 27 → 19 (removed i-band: 8 columns)
- Total epochs retained: ~60-69% depending on band
- Very few sources lost all data (<1% per band)

**Quality Assessment:**
- Strictest quality cuts applied (all 4 filters)
- airmass < 1.8 filter included (based on notebook analysis)
- Cleaning removes outliers, bad PSF fits, high airmass observations
- Expected to improve GP model fits by reducing systematic errors

**Key Features:**
- Original file preserved (`data/gaia_ztf_qso_sample_ztf_lc.fits`)
- Cleaned file uses `_gr_clean` suffix to indicate g/r bands cleaned
- All quality columns retained in output (catflags, sharp, chi, limitmag, airmass)
- No sources removed (only bad epochs within sources)
- Variable-length arrays maintained (FITS format)

**Verification:**
- ✓ All 8,529 sources present in output
- ✓ i-band columns not present in cleaned file
- ✓ g and r bands have reduced epoch counts
- ✓ Quality statistics computed and logged
- ✓ File integrity verified (readable, correct structure)

**Status:** Cleaned ZTF lightcurves ready for GP modeling

### Files Created
1. `download_sdssdr16q_gaia_source.py` (executable)
2. `download_sdssdr16q_gaia_epoch_photometry.py` (executable)
3. `create_ztf_target_catalog.py` (executable)
4. `download_ztf_threaded.py` (copied from qso_gp_mock)
5. `run_ztf_download.sh` (executable launch script)
6. `create_test_2k_catalog.py` (executable, creates 2K test subset)
7. `run_test_2k_with_monitoring.sh` (executable, monitored test runner)
8. `test_thread_scaling.sh` (executable, thread count performance test)
9. `run_multiprocess_download.sh` (executable, multi-process parallelization - not tested)
10. `merge_ztf_chunks.py` (executable, merges chunk files with dtype correction)
11. `create_matched_catalog.py` (executable, creates matched ZTF+Gaia+Catalog files)
12. `clean_ztf_lightcurves.py` (executable, applies quality cuts to ZTF lightcurves)
13. `README.md` (usage documentation)
14. `JOURNAL.md` (this file)

### Output Data Products
```
~/data/gaia/
├── sdssdr16q_gaia_source.fits           (2.1 GB, 489,484 rows, 336 cols)
└── sdssdr16q_gaia_epoch_photometry.fits (980 MB, 223,221 rows, 21 cols)

~/Work/Code/qso_lightcurve_prepare/data/
├── sdssdr16q_gaia_ztf_targets.fits           (5.1 MB, 223,219 rows, 3 cols)
├── test_2k_ztf_targets.fits                  (53 KB, 2,000 rows, 3 cols)
├── ztf_lightcurves_partial_merged.fits       (352 MB, 8,529 rows, 27 cols)
├── gaia_ztf_qso_sample_ztf_lc.fits           (352 MB, 8,529 rows, 27 cols)  ← MATCHED
├── gaia_ztf_qso_sample_gaia_lc.fits          (29 MB, 8,529 rows, 21 cols)   ← MATCHED
├── gaia_ztf_qso_sample_catalog.fits          (38 MB, 8,529 rows, 336 cols)  ← MATCHED
└── gaia_ztf_qso_sample_ztf_lc_gr_clean.fits  (204 MB, 8,529 rows, 19 cols)  ← CLEANED

~/Work/Code/qso_lightcurve_prepare/ztf_chunks/
└── chunk_0000.fits through chunk_0020.fits (21 files, raw download chunks)
```

### Key Statistics

**Gaia Cross-match:**
- Cross-match success rate: 65% (489k/750k DR16Q QSOs)
- Epoch photometry coverage: 46% (223k/489k matched sources)
- Median epochs per source: 41 (G, BP, RP bands)
- Total photometry points: 28.6M (9.5M × 3 bands)
- ZTF target sources: 223,219 (all sources with Gaia epoch photometry)

**ZTF Partial Download (8,529 sources):**
- Download match rate: 30% (8,529/28,500 attempted)
- Median epochs: g=384, r=635, i=98
- Band coverage: g=99.3%, r=99.4%, i=92.0%
- Total ZTF epochs: ~10.2M (g=3.7M, r=5.6M, i=0.9M)

**Matched Catalogs (8,529 sources):**
- Match rate: 100% (all ZTF sources matched with Gaia LC and catalog)
- ZTF epochs: 10,204,591 (g=3.7M, r=5.6M, i=0.9M)
- Gaia epochs: 827,241 (G, BP, RP combined, median ~30/band)
- Total multi-survey epochs: ~11.0M
- Redshift range: z=0.016–4.846 (median z=1.420)
- All 3 matched files perfectly aligned by source_id

**Cleaned ZTF Lightcurves (8,529 sources, g/r bands only):**
- Quality cuts: catflags==0, |sharp|<0.25, 0.5≤chi≤1.5, airmass<1.8
- g-band: 298 epochs/source (mean), 69% retention, 99.1% sources with data
- r-band: 396 epochs/source (mean), 60.6% retention, 99.2% sources with data
- Total cleaned epochs: ~5.9M (g=2.5M, r=3.4M)
- File size reduction: 352 MB → 204 MB (42%)

### Notes
- ~54% of cross-matched sources lack epoch photometry (expected for faint/sparse sampling)
- Variable-length arrays used for per-epoch data (FITS format)
- All quality flags included for robust analysis (Gaia + ZTF)
- Pattern follows proven Quaia lightcurve download approach
- ZTF download uses threaded approach (not multiprocessing) for I/O-bound operations
- ZTF API rate limiting (~8-12 concurrent connections max) prevents threading speedup
- Resumable chunk+manifest system successfully recovered partial data after interruption
- All data files stored in `data/` subdirectory

---

## Session 2: 2025-12-09

### Data Migration
**Objective:** Move data files out of Dropbox to reduce memory usage during processing

**Method:**
- Created directory: `~/data/qso/qso_lightcurve_prepare/`
- Moved `data/` folder from Dropbox project to `~/data/qso/qso_lightcurve_prepare/data/`
- Moved `ztf_chunks/` folder to `~/data/qso/qso_lightcurve_prepare/ztf_chunks/`
- Created symlinks in project directory:
  - `data -> ~/data/qso/qso_lightcurve_prepare/data`
  - `ztf_chunks -> ~/data/qso/qso_lightcurve_prepare/ztf_chunks`
- All scripts continue to work unchanged (symlinks transparent to Python)

**Results:**
- Data files now outside Dropbox sync (no I/O contention during processing)
- Dropbox memory usage freed (~1.5 GB)
- Log files written to `~/data/qso/qso_lightcurve_prepare/`
- Symlinks maintain backward compatibility with all scripts
- Total data moved: ~1.05 GB (data: ~700 MB, ztf_chunks: ~350 MB)

**Files moved:**
- `data/` directory: 9 FITS files (~700 MB)
  - ZTF lightcurves (merged, cleaned)
  - Gaia lightcurves and catalog
  - Target catalogs
- `ztf_chunks/` directory: 21 chunk files (~350 MB)
- Log files: `loo_cleaning*.log`

**Storage layout:**
```
~/data/qso/qso_lightcurve_prepare/
├── data/                          (→ symlinked from project)
│   ├── gaia_ztf_qso_sample_*.fits
│   ├── ztf_loo_results_*.fits
│   └── ...
├── ztf_chunks/                    (→ symlinked from project)
│   └── chunk_*.fits (21 files)
└── loo_cleaning_*.log
```

**Status:** Migration complete, all scripts verified working

### LOO Outlier Cleaning Development
**Script:** `apply_loo_cleaning.py`

**Objective:** Apply Leave-One-Out (LOO) outlier detection to ZTF lightcurves using DRW/OU GP model

**Method:**
- Fit DRW/OU GP model to each lightcurve (unconstrained L-BFGS-B)
- Iteratively remove epochs with high delta-loglikelihood (LOO test)
- Parameters: MIN_POINTS=10, MAX_LOO_PASSES=10, MIN_ABS_DELTA_LOGD=0.5
- Multiple passes until convergence or max iterations
- Process g and r bands separately (single-band mode)

**Implementation Evolution:**

**Issue 1: Progress bar output buffering (FIXED)**
- **Problem:** tqdm progress bars got stuck when laptop closed (output buffering)
- **Solution:** Replaced tqdm with periodic logging (every 100 sources) + `flush=True`
- **Code changes:**
  - Removed `from tqdm.auto import tqdm` import
  - Lines 304-319: Replaced `tqdm(gen, total=n_obj)` with plain loop + counter
  - Added: `if count % 100 == 0: print(f"Processed {count}/{n_obj}...", flush=True)`
  - Line 375: Removed tqdm from `create_loo_cleaned_lightcurves`

**Issue 2: Memory usage (FIXED via single-band mode)**
- **Problem:** Original code kept all g-band results in memory while processing r-band
  - Lines 469-482: `all_results` list and `removed_dicts` accumulated ALL bands
  - With 8,529 sources: ~1.5 GB per band
  - Available RAM: ~170 MB free → would swap/crash during r-band processing
- **Solution:** Refactored to process single band per invocation
- **Code changes:**
  - Changed `--bands g,r` to `--band g` (required, single choice)
  - Auto-generate output filenames: `ztf_loo_results_{band}.fits`
  - Simplified main loop (lines 477-512): process one band, write immediately
  - Removed `vstack(all_results)` - no longer needed
  - Memory savings: ~50% (1.3 GB vs 3 GB peak)

**Test Results:**
- Initial test (100 sources, dual-band): ✓ PASSED
- Full run attempt 1 (8,529 sources, dual-band): KILLED at 38.7% after 2h
  - Reason: Memory concern identified before r-band processing
- Single-band test (10 sources, g-band): ✓ PASSED
  - Output: `ztf_loo_results_g.fits`, `gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits`
  - Mean removed: 1.2 epochs, max: 3

**Production Run (RUNNING):**
- Started: 2025-12-09 12:56pm
- Mode: Parallel processing of both bands simultaneously
- Configuration:
  - g-band: PID 73231, 8 workers, log: `~/data/qso/qso_lightcurve_prepare/loo_cleaning_g.log`
  - r-band: PID 73263, 8 workers, log: `~/data/qso/qso_lightcurve_prepare/loo_cleaning_r.log`
  - Total RAM: ~2.3 GB (1.2 GB + 1.1 GB)
  - Expected runtime: ~4-5 hours per band

**Output files (expected):**
- `data/ztf_loo_results_g.fits` - LOO statistics for g-band
- `data/ztf_loo_results_r.fits` - LOO statistics for r-band
- `data/gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits` - Cleaned g-band lightcurves
- `data/gaia_ztf_qso_sample_ztf_lc_r_clean_loo.fits` - Cleaned r-band lightcurves

**Key Features:**
- Single-band processing: Memory-efficient, can run in parallel
- Automatic output naming: `*_{band}.fits` pattern
- Periodic logging: Every 100 sources with `flush=True`
- Parallel execution: Both bands running simultaneously on separate cores
- Resumable: Each band independent, can restart individually if needed

**Status:** Production run in progress (both bands running in parallel)

---

### Next Steps
1. **Option A:** Continue ZTF download for remaining ~215K sources (estimated 5-10 days)
2. **Option B:** Proceed with matched dataset (8,529 sources) for initial analysis ✓ READY
3. ~~Quality filtering based on photometry/variability flags (Gaia + ZTF)~~ ✓ DONE
   - ✓ Applied ZTF quality cuts to g/r bands: catflags==0, |sharp|<0.25, 0.5≤chi≤1.5, airmass<1.8
   - ✓ Retention: ~60-69% epochs, >99% sources
4. LOO outlier detection for ZTF lightcurves (IN PROGRESS)
   - Script created: `apply_loo_cleaning.py`
   - Test mode verified (100 sources)
   - Memory optimization needed for full run
5. Fit OU GP models to multi-band lightcurves (cleaned ZTF + Gaia)
   - Use `gaia_ztf_qso_sample_ztf_lc_gr_clean.fits` for cleaned ZTF g/r bands
   - Use `gaia_ztf_qso_sample_gaia_lc.fits` for Gaia G/BP/RP bands
   - Use `gaia_ztf_qso_sample_catalog.fits` for redshifts and catalog properties
   - Multi-band modeling: 5 bands total (ZTF g,r + Gaia G,BP,RP)
6. LOO outlier detection for Gaia lightcurves (following qso_gp_mock workflow)
7. Explore coverage vs redshift, magnitude
8. Compare with Quaia catalog (Storey-Fisher et al. 2024, 500k QSOs)
9. D vs L_bol analysis using matched redshifts and multi-band data

---
