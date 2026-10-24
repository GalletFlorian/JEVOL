"""
Stellar evolution model loader.
Reads stellar model files and computes moments of inertia and breakup velocities.
"""

import os

import constants as c
import numpy as np

NSTEP_MAX = 20000


def load_stellar_model(model_file):
    """
    Load a stellar evolution model file.

    File format (after 1-line header):
        idx, Lum, Teff, tstar(yr), rstar(Rsol), k2conv, k2rad, Mrad, Rrad

    Returns arrays (all length n, where n = number of valid lines):
        tstar: age in seconds
        rstar: stellar radius in cm
        lum: luminosity (solar units)
        teff: effective temperature (K)
        k2conv: convective壳 constant
        k2rad: radiative core constant
        mrad: core mass (solar units, relative)
        rrad: core radius (relative to rstar, then converted to cm)
        k2: total k2 = k2conv + k2rad
    """
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")

    tstar = []
    rstar = []
    lum = []
    teff = []
    k2conv = []
    k2rad = []
    mrad = []
    rrad = []

    with open(model_file) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                vals = [float(x) for x in parts[:9]]
            except ValueError:
                continue
            lum.append(vals[1])
            teff.append(vals[2])
            tstar.append(vals[3])
            rstar.append(vals[4])
            k2conv.append(vals[5])
            k2rad.append(vals[6])
            mrad.append(vals[7])
            rrad.append(vals[8])
            if len(tstar) >= NSTEP_MAX:
                break

    n = len(tstar)
    if n == 0:
        raise ValueError(f"No valid data lines in model file: {model_file}")

    tstar = np.array(tstar)
    rstar = np.array(rstar)
    lum = np.array(lum)
    teff = np.array(teff)
    k2conv = np.array(k2conv)
    k2rad_arr = np.array(k2rad)
    mrad = np.array(mrad)
    rrad = np.array(rrad)

    k2 = k2conv + k2rad_arr

    return n, tstar, rstar, lum, teff, k2conv, k2rad_arr, mrad, rrad, k2


def compute_structure(tstar_yr, rstar_rsol, k2conv, k2rad, mrad, rrad_rel, mstar):
    """
    Convert model quantities to CGS and compute derived structure quantities.

    Parameters are raw arrays from load_stellar_model.
    Returns a dict of numpy arrays.
    """
    len(tstar_yr)

    tstar = tstar_yr * c.YEAR
    rstar = rstar_rsol * c.RSOL
    rrad = rrad_rel * rstar  # relative to rstar

    Ienv = k2conv * mstar * c.MSOL * rstar**2
    Icore = k2rad * mstar * c.MSOL * rstar**2

    wbreakup = ((2.0 / 3.0) ** (3.0 / 2.0)) * (c.G * mstar * c.MSOL) ** 0.5 / (rstar**1.5)
    vbreakup = (2.0 * c.G * mstar * c.MSOL / (3.0 * rstar)) ** 0.5

    return {
        "tstar": tstar,
        "rstar": rstar,
        "rrad": rrad,
        "Ienv": Ienv,
        "Icore": Icore,
        "wbreakup": wbreakup,
        "vbreakup": vbreakup,
        "k2conv": k2conv,
        "k2rad": k2rad,
        "mrad": mrad,
    }


def find_indices(tstar, k2rad, tdisk_sec, mstar, n):
    """
    Find critical evolution phase indices.

    Returns:
        index_disk: last step where tstar <= tdisk
        index_dec: first step where k2rad crosses above 1e-9 (core appears)
        index_ms: last step where tstar <= 6.3e7 yr
        index_sun: last step where tstar <= 4.5e9 yr
        ndebut: starting index for output
    """
    tdisk = tdisk_sec

    index_disk = 0
    for j in range(n):
        if tstar[j] <= tdisk:
            index_disk = j

    index_dec = n - 1
    for j in range(1, n):
        if k2rad[j] > 1e-9 and k2rad[j - 1] < 1e-9:
            index_dec = j
            break

    index_ms = 0
    for j in range(n):
        if tstar[j] <= 6.3e7 * c.YEAR:
            index_ms = j

    index_sun = 0
    for j in range(n):
        if tstar[j] <= c.AGESUN * c.YEAR:
            index_sun = j

    if mstar in (1.0, 0.5):
        ndebut = 0
        for j in range(n):
            if tstar[j] < 4.8e5 * c.YEAR:
                ndebut = j
    else:
        ndebut = 0
        for j in range(n):
            if tstar[j] <= 3.0e5 * c.YEAR:
                ndebut = j

    if index_dec == n - 1 and k2rad[n - 1] < 1e-9:
        index_dec = n

    return index_disk, index_dec, index_ms, index_sun, ndebut
