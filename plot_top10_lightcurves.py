#!/usr/bin/env python
"""
Plot the top N SDSS QSOs with the most ZTF data points.
Shows ZTF g, ZTF r, and Gaia G bands on the same plot for each QSO.
"""

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import os

# Set plot style
plt.style.use('default')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

def plot_qso_lightcurves(qso_name, ztf_g_row, ztf_r_row, gaia_row, cat_row, output_dir):
    """Plot all three bands for one QSO."""

    fig, ax = plt.subplots(1, 1, figsize=(14, 6))

    # Get redshift for rest-frame conversion
    z = cat_row['z']

    # Plot ZTF g-band
    if ztf_g_row is not None and len(ztf_g_row['mjd']) > 0:
        mjd_g = ztf_g_row['mjd']
        mag_g = ztf_g_row['mag']
        magerr_g = ztf_g_row['magerr']
        ax.errorbar(mjd_g, mag_g, yerr=magerr_g, fmt='o', markersize=3,
                   alpha=0.6, label=f'ZTF g ({len(mjd_g)} pts)', color='green')

    # Plot ZTF r-band
    if ztf_r_row is not None and len(ztf_r_row['mjd']) > 0:
        mjd_r = ztf_r_row['mjd']
        mag_r = ztf_r_row['mag']
        magerr_r = ztf_r_row['magerr']
        ax.errorbar(mjd_r, mag_r, yerr=magerr_r, fmt='s', markersize=3,
                   alpha=0.6, label=f'ZTF r ({len(mjd_r)} pts)', color='red')

    # Plot Gaia G-band
    if gaia_row is not None and len(gaia_row['g_transit_time']) > 0:
        # Convert Gaia time to MJD
        # g_transit_time is in JD(TCB) - 2455197.5
        # MJD = JD - 2400000.5, so MJD = (g_transit_time + 2455197.5) - 2400000.5 = g_transit_time + 55197
        mjd_gaia = gaia_row['g_transit_time'] + 55197.0
        mag_gaia = gaia_row['g_transit_mag']
        magerr_gaia = gaia_row['g_transit_flux_error'] / gaia_row['g_transit_flux'] * 1.0857

        # Filter valid points
        valid = np.isfinite(mag_gaia) & np.isfinite(magerr_gaia) & (magerr_gaia > 0)
        mjd_gaia = mjd_gaia[valid]
        mag_gaia = mag_gaia[valid]
        magerr_gaia = magerr_gaia[valid]

        ax.errorbar(mjd_gaia, mag_gaia, yerr=magerr_gaia, fmt='D', markersize=4,
                   alpha=0.6, label=f'Gaia G ({len(mjd_gaia)} pts)', color='blue')

    # Formatting
    ax.invert_yaxis()
    ax.set_xlabel('MJD', fontsize=12)
    ax.set_ylabel('Magnitude', fontsize=12)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Title with QSO info
    title = f'{qso_name}  (z={z:.3f})'
    if 'phot_g_mean_mag' in cat_row.colnames:
        gaia_mean_g = cat_row['phot_g_mean_mag']
        title += f'  Gaia G={gaia_mean_g:.2f}'
    ax.set_title(title, fontsize=14, fontweight='bold')

    plt.tight_layout()

    # Save
    safe_name = qso_name.replace('+', 'p').replace('-', 'm').replace('.', '_')
    output_path = os.path.join(output_dir, f'{safe_name}_lightcurve.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return output_path


def main(n_plots=10):
    print("="*70)
    print(f"Plot Top {n_plots} QSOs with Most ZTF Data Points")
    print("="*70)

    # Create output directory
    output_dir = 'plots_top10_lightcurves'
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    print("\nLoading data...")
    ztf_g = Table.read('data/ztf_g_for_gaia_primary.fits')
    ztf_r = Table.read('data/ztf_r_for_gaia_primary.fits')
    gaia = Table.read('data/gaia_g_for_ztf_primary.fits')
    cat = Table.read('data/sdss_gaia_catalog_for_ztf_gaia_lc.fits')

    print(f"  ZTF g-band: {len(ztf_g):,} sources")
    print(f"  ZTF r-band: {len(ztf_r):,} sources")
    print(f"  Gaia: {len(gaia):,} sources")
    print(f"  Catalog: {len(cat):,} sources")

    # Find QSOs with most total ZTF points (g + r)
    print("\nFinding QSOs with most ZTF data points...")

    # Create lookup dictionaries
    ztf_g_dict = {str(name).strip(): i for i, name in enumerate(ztf_g['sdss_name'])}
    ztf_r_dict = {str(name).strip(): i for i, name in enumerate(ztf_r['sdss_name'])}
    gaia_dict = {str(name).strip(): i for i, name in enumerate(gaia['sdss_name'])}
    cat_dict = {str(name).strip(): i for i, name in enumerate(cat['sdss_name'])}

    # Count total points for each QSO
    qso_counts = []
    for name in set(ztf_g_dict.keys()) | set(ztf_r_dict.keys()):
        n_g = len(ztf_g['mjd'][ztf_g_dict[name]]) if name in ztf_g_dict else 0
        n_r = len(ztf_r['mjd'][ztf_r_dict[name]]) if name in ztf_r_dict else 0
        total = n_g + n_r

        # Only include if we have all three bands
        if name in gaia_dict and name in cat_dict and total > 0:
            qso_counts.append((name, total, n_g, n_r))

    # Sort by total points
    qso_counts.sort(key=lambda x: x[1], reverse=True)

    print(f"  Found {len(qso_counts):,} QSOs with data in all bands")
    print(f"\nTop {n_plots} QSOs by total ZTF points:")
    for i, (name, total, n_g, n_r) in enumerate(qso_counts[:min(n_plots, 20)], 1):
        print(f"  {i:2d}. {name:20s}  Total: {total:4d}  (g: {n_g:4d}, r: {n_r:4d})")
    if n_plots > 20:
        print(f"  ... ({n_plots - 20} more)")

    # Plot top N
    print(f"\nGenerating {n_plots} plots...")
    for i, (qso_name, total, n_g, n_r) in enumerate(qso_counts[:n_plots], 1):
        # Get data for this QSO
        ztf_g_row = ztf_g[ztf_g_dict[qso_name]] if qso_name in ztf_g_dict else None
        ztf_r_row = ztf_r[ztf_r_dict[qso_name]] if qso_name in ztf_r_dict else None
        gaia_row = gaia[gaia_dict[qso_name]]
        cat_row = cat[cat_dict[qso_name]]

        output_path = plot_qso_lightcurves(qso_name, ztf_g_row, ztf_r_row, gaia_row,
                                           cat_row, output_dir)
        if i <= 20 or i % 10 == 0:
            print(f"  {i:3d}/{n_plots} {qso_name:20s} → {os.path.basename(output_path)}")

    print(f"\n{'='*70}")
    print("COMPLETE")
    print(f"{'='*70}")
    print(f"\nPlots saved to: {output_dir}/")
    print(f"Total plots: {n_plots}")
    print()


if __name__ == '__main__':
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    main(n)
