"""
DRW/OU GP fitting utilities.

Copied from qso_gp_mock/detect_outliers_loo_bfgs.py
"""

import numpy as np
from celerite2 import terms, GaussianProcess
from scipy import optimize


def fit_drw_unconstrained(times_rest, mags_centered, mag_errs):
    """
    Fit DRW/OU GP with unconstrained L-BFGS-B (no priors).

    Parameters
    ----------
    times_rest : array
        Rest-frame times (days)
    mags_centered : array
        Centered magnitudes (median subtracted)
    mag_errs : array
        Magnitude errors

    Returns
    -------
    dict or None
        Dictionary with keys 'log_sigma', 'log_tau', 'sigma', 'tau', 'success'
        or None if fit failed
    """
    try:
        def neg_ll(params):
            log_sigma, log_tau = params
            sigma = np.exp(log_sigma)
            tau = np.exp(log_tau)

            sigma2 = sigma**2
            kernel = terms.RealTerm(a=sigma2, c=1.0 / tau)
            gp = GaussianProcess(kernel)
            try:
                gp.compute(times_rest, yerr=mag_errs)
                ll = gp.log_likelihood(mags_centered)
            except Exception:
                # Any numerical issue in the GP gets penalised
                return np.inf

            # Extra safety: if ll is NaN or +/-inf, penalise it
            if not np.isfinite(ll):
                return np.inf

            return -ll

        # Initial guess
        y_var = np.nanvar(mags_centered)
        log_sigma_init = 0.5 * np.log(max(y_var, 1e-6))
        baseline = np.max(times_rest) - np.min(times_rest)
        log_tau_init = np.log(max(baseline * 0.1, 1.0))

        soln = optimize.minimize(
            neg_ll,
            [log_sigma_init, log_tau_init],
            method="L-BFGS-B",
            bounds=[(np.log(0.001), np.log(100.0)), (np.log(0.1), np.log(1e4))],
        )

        if not soln.success:
            return None

        log_sigma_opt, log_tau_opt = soln.x
        return {
            "log_sigma": log_sigma_opt,
            "log_tau": log_tau_opt,
            "sigma": np.exp(log_sigma_opt),
            "tau": np.exp(log_tau_opt),
            "success": True,
        }
    except Exception:
        return None
