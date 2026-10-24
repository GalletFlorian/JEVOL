"""
Magnetic braking calculation.
"""

import math

import constants as c
from mass_loss import mdot_calculation


def magnetic_braking(w0, r, lum, t1, t2, ic, config, teff):
    """
    Calculate angular velocity change due to magnetic braking.

    Parameters:
        w0: current angular velocity (s^-1)
        r: stellar radius (cm)
        lum: luminosity (solar units)
        t1: current time (s)
        t2: previous time (s)
        ic: envelope moment of inertia (g cm^2)
        config: simulation config dict
        teff: effective temperature (K)

    Returns:
        Wdownconv: angular velocity lost to braking (s^-1)
    """
    dt = t1 - t2
    mstar = config["mstar"]

    wdown = 0.0

    ksk = config["ksk"]
    kmm = config["kmm"]
    ksc = config["ksc"]
    kmp = config["kmp"]
    brklaw = config["brklaw"]

    wcrit = 0.0
    wsat = 30.0 * c.WSOL

    rsol_ratio = r / c.RSOL
    msol_inv_sqrt = mstar ** (-0.5)

    # Case 1: Skumanich law (dOmega/dt ~ -Omega^3)
    if (wcrit == 0.0 and ksk != 0.0) or (wcrit != 0.0 and w0 < wcrit):
        wdown = ksk * w0**3 * rsol_ratio**0.5 * msol_inv_sqrt * dt / ic

    # Case 2: Mayor-Mermilliod law (dV/dt ~ -V^2)
    elif (wcrit == 0.0 and kmm != 0.0) or (wcrit != 0.0 and w0 >= wcrit):
        if config.get("ksat", 0.0) == 0.0 or w0 < wsat:
            wdown = kmm * w0**2 * rsol_ratio**0.5 * msol_inv_sqrt * dt / ic
        else:
            ksat = config.get("ksat", 0.0)
            wdown = ksat * w0 * rsol_ratio**0.5 * msol_inv_sqrt * dt / ic

    # Case 3: Schatzman law (dV/dt ~ -V^(7/3))
    elif wcrit == 0.0 and ksc != 0.0:
        wdown = ksc * w0 ** (7.0 / 3.0) * rsol_ratio**0.5 * msol_inv_sqrt * dt / ic

    # Case 4: Matt et al. 2012 / Reville / Matt 2015 (kmp != 0)
    elif ksk == 0.0 and kmm == 0.0 and kmp != 0.0:
        mdotstar, bstar, ro, fstar = mdot_calculation(w0, r, lum, config, teff)

        if brklaw == 0:
            # Matt et al. 2012
            K1MP = config["K1MP"]
            K2MP = config["K2MP"]
            K = config["K"]
            m = config["m"]
            M_cgs = mstar * c.MSOL

            num = K1MP**2 * kmp * bstar ** (4 * m) * r ** (5 * m + 2) * w0 * mdotstar ** (1 - 2 * m) * dt
            den = ic * (K2MP**2 * 2 * c.G * M_cgs + K * w0**2 * r**3) ** m
            wdown = num / den if den != 0 else 0.0

        elif brklaw in (1, 2):
            # Reville et al. 2015
            K3 = config["K3"]
            K4 = config["K4"]
            mvic = config["mvic"]
            M_cgs = mstar * c.MSOL

            num = (
                K3**2
                * bstar ** (4 * mvic)
                * r ** (5 * mvic + 2)
                * w0
                * mdotstar ** (1 - 2 * mvic)
                * dt
                * (4 * c.PI) ** (4 * mvic)
            )
            den = ic * (2 * c.G * M_cgs + (2 / K4**2) * w0**2 * r**3) ** mvic
            wdown = num / den if den != 0 else 0.0

        elif brklaw == 3:
            # Matt et al. 2015
            tcz = 314.24 * math.exp(-(teff / 1952.5) - (teff / 6250.0) ** 18) + 0.002
            Rosun = 1.96
            chi = 10.0
            Wsun = 2.87e-6
            tczsun = 12.9
            factor = 0.16
            m_param = config["m"]

            ratio = w0 * r**1.5 * (c.G * mstar * c.MSOL) ** (-0.5)
            gamma = (1 + (ratio / 0.072) ** 2) ** 0.5
            T0 = factor * 5.0e31 * (r / c.RSOL) ** 3.1 * mstar**0.5 * gamma ** (2 * m_param)
            p = 2.3

            if ro <= (Rosun / chi):
                wdown = (T0 * chi**p * (w0 / Wsun)) * dt / ic
            else:
                wdown = (T0 * (tcz / tczsun) ** p * (w0 / Wsun) ** (p + 1)) * dt / ic

    # Case 5: Reiners & Mohanty 2012 (all k = 0)
    elif ksk == 0.0 and kmm == 0.0 and kmp == 0.0:
        Wcr = 8.56e-6
        C = 2.66e3
        M_cgs = mstar * c.MSOL
        r_factor = (r**16 / M_cgs**2) ** (1.0 / 3.0)

        if w0 < Wcr:
            wdown = C * ((w0 / Wcr) ** 4 * w0 * r_factor) * dt / ic
        else:
            wdown = C * (w0 * r_factor) * dt / ic

    return wdown
