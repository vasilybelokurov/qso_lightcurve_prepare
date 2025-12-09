# QSO Lightcurve Preparation Pipeline

Tools for downloading SDSS DR16Q × Gaia DR3 cross-matched catalogs and epoch photometry from WSDB.

## Overview

This pipeline creates two FITS files:
1. **Cross-match catalog** (~400 MB): All SDSS DR16Q and Gaia DR3 source columns
2. **Epoch photometry** (~3 GB): Gaia DR3 lightcurves (g, bp, rp bands)

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

- Cross-match catalog: ~400 MB, ~469k rows, ~350 columns
- Epoch photometry: ~3 GB, ~100M photometry points

## Notes

- Step 1 must complete before running Step 2
- Output locations can be configured in each script
- Conflicting column names are prefixed to avoid collisions
- Uses proven `local_join` pattern from Quaia lightcurve downloads
