"""
Sub-stepping interpolation routines for large velocity changes.
"""

import constants as c
from magnetic_braking import magnetic_braking


def _interpolate_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, mstar, j, ninter):
    """
    Linearly interpolate stellar structure between steps j-1 and j.
    """
    ninter1 = ninter + 1
    ti = [0.0] * ninter1
    ri = [0.0] * ninter1
    kici = [0.0] * ninter1
    kiri = [0.0] * ninter1
    Ici = [0.0] * ninter1
    Iri = [0.0] * ninter1
    rri = [0.0] * ninter1
    mri = [0.0] * ninter1

    for k in range(ninter1):
        frac = k / ninter
        ti[k] = tstar[j - 1] + frac * (tstar[j] - tstar[j - 1])
        ri[k] = rstar[j - 1] + frac * (rstar[j] - rstar[j - 1])
        kici[k] = k2conv[j - 1] + frac * (k2conv[j] - k2conv[j - 1])
        kiri[k] = k2rad[j - 1] + frac * (k2rad[j] - k2rad[j - 1])
        rri[k] = rrad[j - 1] + frac * (rrad[j] - rrad[j - 1])
        mri[k] = mrad[j - 1] + frac * (mrad[j] - mrad[j - 1])
        Ici[k] = kici[k] * mstar * c.MSOL * ri[k] ** 2
        Iri[k] = kiri[k] * mstar * c.MSOL * ri[k] ** 2

    return ti, ri, kici, kiri, Ici, Iri, rri, mri


def interpolA(Wdownconv, Wup, W0, Wconv_prev, rstar, lum, teff, tstar, j, mstar, config):
    """
    Sub-step interpolation for Phase II.A (envelope only, no core).

    If the velocity change is >10% of W0, subdivide the timestep into
    ninter sub-steps and integrate forward.

    Returns:
        (Winterfin, Wdownconv_total) - final angular velocity and total braking
    """
    Wtest = max(abs(Wdownconv), abs(Wup))

    if Wtest <= 0.1 * abs(W0):
        return W0 + Wup - Wdownconv, Wdownconv

    ninter = min(max(int(Wtest / (0.1 * abs(W0))), 10), 200)

    k2conv = config.get("_k2conv", None)
    k2rad = config.get("_k2rad", None)
    mrad_arr = config.get("_mrad", None)
    rrad_arr = config.get("_rrad", None)
    rstar_arr = config.get("_rstar", None)
    lum_arr = config.get("_lum", None)
    teff_arr = config.get("_teff", None)
    tstar_arr = config.get("_tstar", None)

    ti, ri, kici, kiri, Ici, Iri, rri, mri = _interpolate_structure(
        tstar_arr, rstar_arr, k2conv, k2rad, mrad_arr, rrad_arr, mstar, j, ninter
    )

    # Interpolate luminosity and Teff between j-1 and j
    Li = [0.0] * (ninter + 1)
    Teffi = [0.0] * (ninter + 1)
    for k in range(ninter + 1):
        frac = k / ninter
        Li[k] = lum_arr[j - 1] + frac * (lum_arr[j] - lum_arr[j - 1])
        Teffi[k] = teff_arr[j - 1] + frac * (teff_arr[j] - teff_arr[j - 1])

    Wi = [0.0] * (ninter + 1)
    Wi[0] = W0
    sum_wdown = 0.0

    for k in range(1, ninter + 1):
        Wup_k = Wi[k - 1] * (Ici[k - 1] / Ici[k] - 1.0) if Ici[k] > 0 else 0.0
        Wdown_k = magnetic_braking(Wi[k - 1], ri[k - 1], Li[k - 1], ti[k], ti[k - 1], Ici[k - 1], config, Teffi[k - 1])
        Wi[k] = Wi[k - 1] + Wup_k - Wdown_k
        sum_wdown += Wdown_k

    return Wi[ninter], sum_wdown


def interpolB(Juprad, dJ, dWr, Wur, Wr0, Wc0, Iri_arr, Ici_arr, rri_arr, mri_arr, ti_arr, j, mstar, config):
    """
    Sub-step interpolation for Phase II.B (core evolving, disk-locked envelope).

    Returns:
        (Winterfinr, taudecfin) - final core angular velocity and coupling timescale
    """
    Wtest = max(abs(Juprad / Iri_arr[j]) if Iri_arr[j] > 0 else 0.0, abs(dWr), abs(Wur))

    if Wtest < 0.1 * abs(Wr0):
        val = Wr0 + Juprad / Iri_arr[j] - dWr + Wur if Iri_arr[j] > 0 else Wr0
        return val, config["taudecinit"]

    ninter = min(max(int(Wtest / (0.1 * abs(Wr0))) if abs(Wr0) > 1e-10 else 10, 10), 200)

    tstar_arr = config.get("_tstar", None)
    rstar_arr = config.get("_rstar", None)
    k2conv = config.get("_k2conv", None)
    k2rad = config.get("_k2rad", None)
    mrad_arr = config.get("_mrad", None)
    rrad_arr = config.get("_rrad", None)

    ti_local, ri_local, kici, kiri, Ici_local, Iri_local, rri_local, mri_local = _interpolate_structure(
        tstar_arr, rstar_arr, k2conv, k2rad, mrad_arr, rrad_arr, mstar, j, ninter
    )

    taudecinit = config["taudecinit"]
    Wri = [0.0] * (ninter + 1)
    Wri[0] = Wr0
    dJ_local = dJ

    for k in range(1, ninter + 1):
        wur_k = Wri[k - 1] * (Iri_local[k - 1] / Iri_local[k] - 1.0) if Iri_local[k] > 0 else 0.0
        Juprad_k = (2.0 / 3.0) * rri_local[k] ** 2 * Wc0 * (mri_local[k] - mri_local[k - 1]) * c.MSOL
        dWr_k = (dJ_local * (ti_local[k] - ti_local[k - 1]) / (taudecinit * Iri_local[k])) if Iri_local[k] > 0 else 0.0
        Wri[k] = Wri[k - 1] + Juprad_k / Iri_local[k] - dWr_k + wur_k if Iri_local[k] > 0 else 0.0

        if Wc0 < Wri[k]:
            dJ_local = Iri_local[k] * Ici_local[k] / (Iri_local[k] + Ici_local[k]) * (Wri[k] - Wc0)
        elif Wc0 == Wri[k]:
            dJ_local = 0.0
        else:
            dJ_local = -(Iri_local[k] * Ici_local[k] / (Iri_local[k] + Ici_local[k]) * (Wc0 - Wri[k]))

    return Wri[ninter], taudecinit


def interpol3(Wdownconv, Juprad, Wc0, Wr0, dWr, dWc, dJ, Wuprad, Wupconv, j, mstar, config):
    """
    Sub-step interpolation for Phase III (both core and envelope evolving).

    Returns:
        (Winterfin, Winterfinr) - final envelope and core angular velocities
    """
    tstar_arr = config.get("_tstar", None)
    rstar_arr = config.get("_rstar", None)
    k2conv = config.get("_k2conv", None)
    k2rad = config.get("_k2rad", None)
    mrad_arr = config.get("_mrad", None)
    rrad_arr = config.get("_rrad", None)
    lum_arr = config.get("_lum", None)
    teff_arr = config.get("_teff", None)
    taudecinit = config["taudecinit"]

    # Always sub-step (Fortran always overrides ninter to 2000)
    ninter = 20

    ti, ri, kici, kiri, Ici, Iri, rri, mri = _interpolate_structure(
        tstar_arr, rstar_arr, k2conv, k2rad, mrad_arr, rrad_arr, mstar, j, ninter
    )

    # Interpolate luminosity and Teff between j-1 and j
    Li = [0.0] * (ninter + 1)
    Teffi = [0.0] * (ninter + 1)
    for k in range(ninter + 1):
        frac = k / ninter
        Li[k] = lum_arr[j - 1] + frac * (lum_arr[j] - lum_arr[j - 1])
        Teffi[k] = teff_arr[j - 1] + frac * (teff_arr[j] - teff_arr[j - 1])

    Wci = [0.0] * (ninter + 1)
    Wri = [0.0] * (ninter + 1)
    Wci[0] = Wc0
    Wri[0] = Wr0
    dJ_local = dJ

    for k in range(1, ninter + 1):
        wuc = Wci[k - 1] * (Ici[k - 1] / Ici[k] - 1.0) if Ici[k] > 0 else 0.0
        wur = Wri[k - 1] * (Iri[k - 1] / Iri[k] - 1.0) if Iri[k] > 0 else 0.0

        Wdown_k = magnetic_braking(Wci[k - 1], ri[k - 1], Li[k - 1], ti[k], ti[k - 1], Ici[k - 1], config, Teffi[k - 1])

        Juprad_k = (2.0 / 3.0) * rri[k] ** 2 * Wci[k - 1] * (mri[k] - mri[k - 1]) * c.MSOL
        dWc_k = (dJ_local * (ti[k] - ti[k - 1]) / (taudecinit * Ici[k])) if Ici[k] > 0 else 0.0
        dWr_k = (dJ_local * (ti[k] - ti[k - 1]) / (taudecinit * Iri[k])) if Iri[k] > 0 else 0.0

        Wci[k] = (Wci[k - 1] - Wdown_k + wuc - Juprad_k / Ici[k] + dWc_k) if Ici[k] > 0 else Wci[k - 1]

        Juprad_k2 = (2.0 / 3.0) * rri[k] ** 2 * Wci[k] * (mri[k] - mri[k - 1]) * c.MSOL
        Wri[k] = (Wri[k - 1] + wur + Juprad_k2 / Iri[k] - dWr_k) if Iri[k] > 0 else 0.0

        if Wci[k] < Wri[k]:
            dJ_local = Iri[k] * Ici[k] / (Iri[k] + Ici[k]) * (Wri[k] - Wci[k])
        elif Wci[k] == Wri[k]:
            dJ_local = 0.0
        else:
            dJ_local = -(Iri[k] * Ici[k] / (Iri[k] + Ici[k]) * (Wci[k] - Wri[k]))

    return Wci[ninter], Wri[ninter]
