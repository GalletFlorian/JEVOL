"""
Main rotational evolution loop.
"""

import constants as c
import numpy as np
from interpolation import interpol3, interpolA, interpolB
from magnetic_braking import magnetic_braking


def compute_taudec(wrad, wconv, taudecinit):
    """
    Compute coupling timescale.
    In the Fortran code, taudec is forced constant (taudecinit) despite
    the Spada formula being implemented. We match this behavior.
    """
    return taudecinit


def run_evolution(struct, config, indices):
    """
    Run the 3-phase rotational evolution.

    Parameters:
        struct: dict from stellar_model.compute_structure
        config: dict from parameters.build_config
        indices: tuple (index_disk, index_dec, index_ms, index_sun, ndebut)

    Returns:
        dict with evolution results
    """
    n = len(struct["tstar"])
    tstar = struct["tstar"]
    rstar = struct["rstar"]
    Ienv = struct["Ienv"]
    Icore = struct["Icore"]
    lum = config.get("_lum", np.ones(n))
    teff_arr = config.get("_teff", np.ones(n) * 5778.0)
    k2conv = struct["k2conv"]
    k2rad = struct["k2rad"]
    mrad = struct["mrad"]
    rrad = struct["rrad"]

    mstar = config["mstar"]
    pinit = config["pinit"]
    taudecinit = config["taudecinit"]

    index_disk, index_dec, index_ms, index_sun, ndebut = indices

    winit = c.DPI / (pinit * c.DAY)

    wconv = np.zeros(n)
    wrad = np.zeros(n)
    dj = np.zeros(n)
    amconv = np.zeros(n)
    amrad = np.zeros(n)
    taudec_arr = np.full(n, taudecinit)

    # Phase I: Disk-coupled, no core
    if index_dec >= index_disk:
        index1 = index_disk
        index2 = max(index_dec - 1, 0)
    else:
        index1 = max(index_dec - 1, 0)
        index2 = index_disk

    for j in range(1, index1 + 1):
        wconv[j] = winit
        wrad[j] = 0.0
        dj[j] = 0.0
        amconv[j] = Ienv[j] * wconv[j]
        amrad[j] = Icore[j] * wrad[j]
        taudec_arr[j] = taudecinit

    # Enrich config with arrays needed by interpolation/magnetic_braking
    config["_rstar"] = rstar
    config["_lum"] = lum
    config["_teff"] = teff_arr
    config["_k2conv"] = k2conv
    config["_k2rad"] = k2rad
    config["_mrad"] = mrad
    config["_rrad"] = rrad
    config["_tstar"] = tstar

    if index_dec >= index_disk:
        # Case A: disk disappears before decoupling
        for j in range(index1 + 1, min(index2 + 1, n)):
            if Ienv[j] > 0 and Ienv[j - 1] > 0:
                wupconv = wconv[j - 1] * (Ienv[j - 1] / Ienv[j] - 1.0)
            else:
                wupconv = 0.0

            wdownconv = magnetic_braking(
                wconv[j - 1], rstar[j - 1], lum[j - 1], tstar[j], tstar[j - 1], Ienv[j - 1], config, teff_arr[j - 1]
            )

            # Sub-step if velocity change >10% of current velocity
            wconv[j], _ = interpolA(
                wdownconv, wupconv, wconv[j - 1], wconv[j - 1], rstar, lum, teff_arr, tstar, j, mstar, config
            )
            wrad[j] = 0.0
            dj[j] = 0.0
            amconv[j] = Ienv[j] * wconv[j]
            amrad[j] = Icore[j] * wrad[j]
            taudec_arr[j] = taudecinit

    else:
        # Case B: disk disappears after decoupling
        deltaJ = 0.0
        for j in range(index1 + 1, min(index2 + 1, n)):
            if mrad[j] / mstar > 0.98:
                config["K1MP"] -= config["K1MP"] * 0.2
            wconv[j] = winit

            if Icore[j] > 0 and Icore[j - 1] > 0:
                wuprad = wrad[j - 1] * (Icore[j - 1] / Icore[j] - 1.0)
            else:
                wuprad = 0.0

            Juprad = (2.0 / 3.0) * rrad[j] ** 2 * wconv[j - 1] * (mrad[j] - mrad[j - 1]) * c.MSOL

            taudec_arr[j] = compute_taudec(wrad[j - 1], wconv[j - 1], taudecinit)

            if Icore[j] > 0:
                deltaWr = deltaJ * (tstar[j] - tstar[j - 1]) / (taudec_arr[j] * Icore[j])
            else:
                deltaWr = 0.0

            # Sub-step if velocity change >10% of current core velocity
            if j > index1 + 1 and Icore[j] > 0:
                wrad[j], _ = interpolB(
                    Juprad,
                    deltaJ,
                    deltaWr,
                    wuprad,
                    wrad[j - 1],
                    wconv[j - 1],
                    Icore,
                    Ienv,
                    rrad,
                    mrad,
                    tstar,
                    j,
                    mstar,
                    config,
                )
            else:
                if Icore[j] > 0:
                    wrad[j] = wrad[j - 1] + wuprad + Juprad / Icore[j] - deltaWr
                else:
                    wrad[j] = 0.0

            dj[j] = deltaJ
            amrad[j] = Icore[j] * wrad[j]
            amconv[j] = Ienv[j] * wconv[j]

            if j == index1 + 1:
                wrad[j] = wconv[j]

            # Compute deltaJ for next step (Fortran lines 664-682)
            if Icore[j] > 0 and Ienv[j] > 0:
                if wconv[j] < wrad[j]:
                    deltaJ = Icore[j] * Ienv[j] / (Icore[j] + Ienv[j]) * (wrad[j] - wconv[j])
                elif wconv[j] == wrad[j]:
                    deltaJ = 0.0
                else:
                    deltaJ = -(Icore[j] * Ienv[j] / (Icore[j] + Ienv[j]) * (wconv[j] - wrad[j]))
            else:
                deltaJ = 0.0

    # Phase III: Core-envelope decoupled, star-disk uncoupled
    for j in range(index2 + 1, n):
        if Ienv[j] > 0 and Ienv[j - 1] > 0:
            wupconv = wconv[j - 1] * (Ienv[j - 1] / Ienv[j] - 1.0)
        else:
            wupconv = 0.0

        if Icore[j] > 0 and Icore[j - 1] > 0:
            wuprad = wrad[j - 1] * (Icore[j - 1] / Icore[j] - 1.0)
        else:
            wuprad = 0.0

        wdownconv = magnetic_braking(
            wconv[j - 1], rstar[j - 1], lum[j - 1], tstar[j], tstar[j - 1], Ienv[j - 1], config, teff_arr[j - 1]
        )

        taudec_arr[j] = compute_taudec(wrad[j - 1], wconv[j - 1], taudecinit)

        Juprad = (2.0 / 3.0) * rrad[j] ** 2 * wconv[j - 1] * (mrad[j] - mrad[j - 1]) * c.MSOL

        if Icore[j] > 0:
            deltaWr = dj[j - 1] * (tstar[j] - tstar[j - 1]) / (taudec_arr[j] * Icore[j])
        else:
            deltaWr = 0.0

        deltaWc = (dj[j - 1] * (tstar[j] - tstar[j - 1]) / (taudec_arr[j] * Ienv[j])) if Ienv[j] > 0 else 0.0

        # Sub-step if velocity change >5% of current velocity
        if j > index_dec:
            wconv[j], wrad[j] = interpol3(
                wdownconv,
                Juprad,
                wconv[j - 1],
                wrad[j - 1],
                deltaWr,
                deltaWc,
                dj[j - 1],
                wuprad,
                wupconv,
                j,
                mstar,
                config,
            )
        else:
            # Envelope update
            if Ienv[j] > 0:
                wconv[j] = wconv[j - 1] + wupconv - wdownconv - Juprad / Ienv[j] + deltaWc
            else:
                wconv[j] = wconv[j - 1]

            # Recompute Juprad with updated wconv
            Juprad = (2.0 / 3.0) * rrad[j] ** 2 * wconv[j] * (mrad[j] - mrad[j - 1]) * c.MSOL

            # Core update
            if Icore[j] > 0:
                wrad[j] = wrad[j - 1] + wuprad + Juprad / Icore[j] - deltaWr
            else:
                wrad[j] = 0.0

        # Force equal rotation at decoupling step
        if j == index_dec:
            wrad[j] = wconv[j]
            taudec_arr[j] = taudecinit

        # Core-envelope angular momentum exchange
        if wconv[j] < wrad[j]:
            dj[j] = Icore[j] * Ienv[j] / (Icore[j] + Ienv[j]) * (wrad[j] - wconv[j])
        elif wconv[j] == wrad[j]:
            dj[j] = 0.0
        else:
            dj[j] = -(Icore[j] * Ienv[j] / (Icore[j] + Ienv[j]) * (wconv[j] - wrad[j]))

        amconv[j] = Ienv[j] * wconv[j]
        amrad[j] = Icore[j] * wrad[j]

    return {
        "wconv": wconv,
        "wrad": wrad,
        "dj": dj,
        "amconv": amconv,
        "amrad": amrad,
        "taudec": taudec_arr,
        "winit": winit,
    }
