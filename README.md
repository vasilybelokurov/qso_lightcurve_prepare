# QSO Lightcurve Preparation Pipeline

Tools for downloading SDSS DR16Q × Gaia DR3 cross-matched catalogs and epoch photometry from WSDB, and preparing quality-controlled lightcurves for variability analysis.

## Quick Reference: Data Files Summary

| Dataset | Location | Size | N Sources | Bands | Status | Description |
|---------|----------|------|-----------|-------|--------|-------------|
| **GAIA DR3 × SDSS DR16Q** |
| Gaia source catalog | `~/data/gaia/sdssdr16q_gaia_source.fits` | 2.1 GB | 489,484 | - | Raw | Full cross-match with all Gaia+SDSS columns |
| Gaia epoch photometry | `~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits` | 980 MB | 223,221 | G,BP,RP | Raw | Epoch photometry (2014-07-25 to 2017-05-28, 34 months) |
| **ZTF DR19 × SDSS DR16Q** |
| ZTF primary g-band | `~/data/ztf/qso/ztf_g_band_lc_primary.h5` | 3.6 GB | 706,985 | g | Deduplicated + catflags=0 | Primary lightcurves (2018-03 to 2024-10, 6.6 years) |
| ZTF primary r-band | `~/data/ztf/qso/ztf_r_band_lc_primary.h5` | 4.9 GB | 722,701 | r | Deduplicated + catflags=0 | Primary lightcurves (2018-03 to 2024-10, 6.6 years) |
| ZTF primary catalog | `data/ztf_primary_sources_catalog.fits` | 56 MB | 731,702 | - | - | Combined g+r sdss_name list |
| **MATCHED ZTF+GAIA DATASET (220k sources)** |
| SDSS+Gaia catalog | `data/sdss_gaia_catalog_for_ztf_gaia_lc.fits` | 970 MB | 219,830 | - | Raw | Full SDSS+Gaia properties (336 cols) |
| SDSS+Gaia+Wushen VAC 2024 catalog | `data/sdss_gaia_catalog_for_ztf_gaia_lc_with_wushen_vac_2024.fits` | 1.8 GB | 219,830 | - | Augmented | Adds `wushen_*` columns from `sdssdr16qso.wushen_vac_2024` |
| Gaia lightcurves | `data/gaia_g_for_ztf_primary.fits` | 974 MB | 219,778 | G,BP,RP | Raw | Gaia epoch photometry |
| ZTF g-band lightcurves | `data/ztf_g_for_gaia_primary.fits` | 1.7 GB | 217,615 | g | Deduplicated + catflags=0 | ZTF g-band epoch photometry |
| ZTF r-band lightcurves | `data/ztf_r_for_gaia_primary.fits` | 2.3 GB | 218,742 | r | Deduplicated + catflags=0 | ZTF r-band epoch photometry |
| **ALL ZTF PRIMARY SOURCES (482k sources)** |
| SDSS+Gaia catalog | `data/sdss_gaia_catalog_for_ztf_primary.fits` | 2.1 GB | 482,161 | - | Raw | Full SDSS+Gaia properties (336 cols) |
| **LOCAL WORKING SUBSET (8.5k sample)** |
| Catalog | `data/gaia_ztf_qso_sample_catalog.fits` | 38 MB | 8,529 | - | - | Properties for test sample |
| Gaia G LOO cleaned | `data/gaia_ztf_qso_sample_gaia_lc_g_clean_loo.fits` | 29 MB | 8,529 | G | LOO cleaned | Test sample |
| Gaia BP LOO cleaned | `data/gaia_ztf_qso_sample_gaia_lc_bp_clean_loo.fits` | 28 MB | 8,488 | BP | LOO cleaned | Test sample |
| Gaia RP LOO cleaned | `data/gaia_ztf_qso_sample_gaia_lc_rp_clean_loo.fits` | 28 MB | 8,502 | RP | LOO cleaned | Test sample |
| ZTF g LOO cleaned | `data/gaia_ztf_qso_sample_ztf_lc_g_clean_loo.fits` | 88 MB | 8,371 | g | LOO cleaned | Test sample |
| ZTF r LOO cleaned | `data/gaia_ztf_qso_sample_ztf_lc_r_clean_loo.fits` | 116 MB | 8,449 | r | LOO cleaned | Test sample |

**Key identifiers:**
- Primary key: `sdss_name` (format: `HHMMSS.ss+DDMMSS.s`)
- Gaia: `source_id` (int64)
- ZTF: `objectid` (int64)

**Time formats:**
- **Gaia times** (`g_transit_time`, `bp_obs_time`, `rp_obs_time`): JD(TCB) - 2455197.5 days
  - Reference epoch: JD 2455197.5 (TCB) = 2010-01-01T00:00:00
  - **Conversion to MJD**: `MJD = g_transit_time + 55197`
  - Documentation: [Gaia DR3 epoch_photometry](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_photometry/ssec_dm_epoch_photometry.html)
- **ZTF times** (`mjd`): Standard MJD (UTC)

**Storage locations:**
- Full datasets: `~/data/gaia/` and `~/data/ztf/qso/`
- Working files: `data/` (symlinked to `~/data/qso/qso_lightcurve_prepare/data/`)

**Time coverage:**
- **Gaia DR3:** 2014-07-25 to 2017-05-28 (34 months) [[Gaia SSDC](https://gaia.ssdc.asi.it/news.php), [Eyer et al. 2023](https://www.aanda.org/articles/aa/full_html/2023/06/aa44242-22/aa44242-22.html)]
- **ZTF DR19:** 2018-03-17 to 2024-10-31 (6.6 years) [[ZTF Public Releases](https://www.ztf.caltech.edu/ztf-public-releases.html), [IRSA](https://irsa.ipac.caltech.edu/Missions/ztf.html)]
- **Note:** No temporal overlap between Gaia DR3 and ZTF DR19 in this dataset

## Overview

This pipeline provides:
1. **Cross-match catalogs**: SDSS DR16Q × Gaia DR3 and SDSS DR16Q × ZTF DR19
2. **Epoch photometry**: Gaia DR3 (G, BP, RP) and ZTF (g, r) lightcurves
3. **Quality control**: LOO outlier detection and cleaning

## Requirements

- Python 3.x
- `sqlutilpy` (for WSDB access)
- `astropy`
- WSDB credentials configured

## Usage

### Step 1: Download Cross-Match Catalog

```bash
python download_sdssdr16q_gaia_source.py
```

**Output:** Default location (configure in script)

**What it does:**
- Spatially cross-matches SDSS DR16Q final QSOs with Gaia DR3 sources (1 arcsec radius)
- Downloads **all columns** from both catalogs
- Renames conflicting columns (`ra`, `dec`, `duplicated_source`) with prefixes:
  - SDSS: `sdss_ra`, `sdss_dec`, `sdss_duplicated_source`
  - Gaia: `gaia_ra`, `gaia_dec`, `gaia_duplicated_source`
- Expected: ~469k QSOs matched with Gaia

**Runtime:** ~3-5 minutes

### Step 2: Download Epoch Photometry

```bash
python download_sdssdr16q_gaia_epoch_photometry.py
```

**Output:** Default location (configure in script)

**What it does:**
- Loads source_ids from Step 1 catalog
- Uses `local_join` to upload source_ids to WSDB
- Downloads epoch photometry for all 3 bands (g, bp, rp)
- Includes variability and photometry quality flags

**Runtime:** ~10-30 minutes

## Output Files

**IMPORTANT: Full SDSS DR16Q × Gaia DR3 dataset location:**
- Cross-match catalog: `~/data/gaia/sdssdr16q_gaia_source.fits` (2.1 GB, ~469k sources)
- Epoch photometry: `~/data/gaia/sdssdr16q_gaia_epoch_photometry.fits` (980 MB, 223,221 sources with lightcurves)

**Local working subset (8.5k sources with ZTF data):**
- Catalog: `data/gaia_ztf_qso_sample_catalog.fits` (38 MB, 8,529 sources)
- G-band LOO cleaned: `data/gaia_ztf_qso_sample_gaia_lc_g_clean_loo.fits` (29 MB)
- BP-band LOO cleaned: `data/gaia_ztf_qso_sample_gaia_lc_bp_clean_loo.fits` (28 MB)
- RP-band LOO cleaned: `data/gaia_ztf_qso_sample_gaia_lc_rp_clean_loo.fits` (28 MB)

## Matched ZTF+Gaia Dataset (220k sources)

The **primary analysis-ready dataset** contains ~220k SDSS QSOs with BOTH ZTF and Gaia epoch photometry:

### Files in `data/`:

1. **Catalog:** `sdss_gaia_catalog_for_ztf_gaia_lc.fits` (970 MB, 219,830 sources)
   - All 336 columns from SDSS DR16Q + Gaia DR3 source catalogs
   - Conflicting columns prefixed: `sdss_ra`, `sdss_dec`, `gaia_ra`, `gaia_dec`, etc.
   - Key columns: `sdss_name`, `z`, `source_id`, `phot_g_mean_mag`, `parallax`, `pm`, etc.

   **Optional augmented version (adds Wushen VAC 2024 columns):**
   - Output: `sdss_gaia_catalog_for_ztf_gaia_lc_with_wushen_vac_2024.fits` (219,830 sources, +114 cols)
   - Join key: `sdss_name`
   - Command:
     ```bash
     source ~/Work/venvs/.venv/bin/activate
     python add_wushen_vac_2024_to_catalog.py \
       --input data/sdss_gaia_catalog_for_ztf_gaia_lc.fits \
       --output data/sdss_gaia_catalog_for_ztf_gaia_lc_with_wushen_vac_2024.fits \
       --input-key sdss_name --vac-key sdss_name --prefix-all
     ```

2. **Gaia lightcurves:** `gaia_g_for_ztf_primary.fits` (974 MB, 219,778 sources)
   - G, BP, RP band epoch photometry
   - Variable-length arrays: `g_transit_time`, `g_transit_mag`, `g_transit_flux`, etc.
   - RAW (not LOO-cleaned)

3. **ZTF g-band:** `ztf_g_for_gaia_primary.fits` (1.7 GB, 217,615 sources)
   - Variable-length arrays: `mjd`, `mag`, `magerr`, `clrcoeff`
   - Already filtered: deduplicated (primary objectid) + catflags=0
   - Metadata: `n_epochs_clean`, `n_epochs_raw`, `time_baseline`, `separation_arcsec`

4. **ZTF r-band:** `ztf_r_for_gaia_primary.fits` (2.3 GB, 218,742 sources)
   - Same structure as g-band

**Matching:** All files use `sdss_name` as the common identifier.

**Usage example:**
```python
from astropy.table import Table

# Load catalog
cat = Table.read('data/sdss_gaia_catalog_for_ztf_gaia_lc.fits')

# Load lightcurves
gaia = Table.read('data/gaia_g_for_ztf_primary.fits')
ztf_g = Table.read('data/ztf_g_for_gaia_primary.fits')
ztf_r = Table.read('data/ztf_r_for_gaia_primary.fits')

# Match by sdss_name
idx = (cat['sdss_name'] == gaia['sdss_name'][0])
```

## All ZTF Primary Sources (482k sources)

**Catalog:** `sdss_gaia_catalog_for_ztf_primary.fits` (2.1 GB, 482,161 sources)
- Contains all sources from ZTF primary HDF5 files that have Gaia matches
- Includes sources without Gaia epoch photometry (but with Gaia source properties)
- Use this for analysis requiring the full ZTF sample

## Notes

- Step 1 must complete before running Step 2
- Full datasets stored in `~/data/gaia/` and `~/data/ztf/qso/`
- Working subset stored in `data/` (symlinked to external storage)
- Conflicting column names are prefixed to avoid collisions
- Uses proven `local_join` pattern from Quaia lightcurve downloads
