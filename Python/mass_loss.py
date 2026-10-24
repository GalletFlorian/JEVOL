"""
Mass loss rate and magnetic field calculation.
"""

import math

import constants as c


def mdot_calculation(wstar, rstar, lstar, config, teff):
    """
    Calculate mass loss rate, magnetic field, Rossby number, and filling factor.

    Parameters:
        wstar: stellar angular velocity (s^-1)
        rstar: stellar radius (cm)
        lstar: luminosity (solar units)
        config: simulation config dict
        teff: effective temperature (K)

    Returns:
        (mdotstar, bstar, Ro, fstar) in CGS units
    """
    mstar = config["mstar"]
    M = mstar * c.MSOL
    lstar * c.LSOL

    grav = M * c.G / rstar**2
    gravsol = c.MSOL * c.G / c.RSOL**2
    lgrav = math.log10(grav) if grav > 0 else 0.0

    theta = 1.0 / 3.0
    ratioZ = 1.0
    h = 0.5

    # Rossby number
    period = 2.0 * c.PI / (wstar * c.DAY) if wstar > 0 else 1e30
    tauc = 314.24 * math.exp(-(teff / 1952.5) - (teff / 6250.0) ** 18) + 0.002

    gratio = gravsol / grav if grav > 0 else 1.0
    if gratio > 1.0:
        tauc = tauc * gratio**0.18

    Rosol = 1.96
    Ro = period / tauc if tauc > 0 else 1e30
    x = Ro / Rosol

    # Filling factor (Boussinesq model)
    fmod = 0.55 / (1.0 + (x / 0.16) ** 2.3) ** 1.22
    fstar = fmod
    fstarTR = fstar**theta

    # Density (Cranmer approximation)
    xteff = teff / 1000.0
    C0 = (
        4.0872049
        - 48.979961 * xteff
        + 47.135345 * xteff**2
        - 20.204336 * xteff**3
        + 4.4110828 * xteff**4
        - 0.48112223 * xteff**5
        + 0.020825121 * xteff**6
    )
    C1 = (
        24.562656
        - 25.713078 * xteff
        + 9.8731818 * xteff**2
        - 1.3825986 * xteff**3
        - 0.055190235 * xteff**4
        + 0.031616983 * xteff**5
        - 0.0021585992 * xteff**6
    )

    rhostar = 10.0 ** (C0 + C1 * lgrav) if lgrav > -20 else 1e-30

    # Magnetic field
    kb = 1.380622e-16
    mh = 1.67333e-24
    mu = 7.0 / 4.0 + 0.5 * math.tanh((3500.0 - teff) / 600.0)

    Bstar = 1.13 * math.sqrt(8.0 * c.PI * rhostar * kb * teff / (mu * mh)) if (rhostar > 0 and teff > 0) else 0.0
    Bstar * fstar

    # Alfv\u00e9n velocity
    vesc2 = 2.0 * c.G * M / rstar
    vesc = math.sqrt(vesc2) if vesc2 > 0 else 1.0

    VA = Bstar / math.sqrt(4.0 * c.PI * rhostar) if rhostar > 0 else 1e30

    # Hot wind: energy flux
    F0 = 5.724 * math.exp(-lgrav / 11.48) * 1e9 if lgrav > -20 else 0.0
    T0 = 1000.0 * (5.624 + 0.6002 * lgrav)
    eps = 6.774 + 0.5057 * lgrav
    FA = F0 * (teff / T0) ** eps * math.exp(-((teff / T0) ** 25))
    FA = FA / 2.5  # Best fit for 1 Msol

    lammaxsol = 7.4e-23 + 4.2e-22
    lammax = 1e-23 * (7.4 + 42.0 * ratioZ**1.13)
    lampersol = 3e7
    Hsol = 1.39e7
    lamperstar = (lampersol / Hsol) * kb * teff / (mu * mh * grav) if grav > 0 else 1e30

    vper = math.sqrt(FA / (rhostar * VA)) if (rhostar > 0 and VA > 0) else 1.0
    Qstar = c.ALPHAS if False else 0.5  # alpha parameter
    Qstar = 0.5 * rhostar * vper**3 / lamperstar if lamperstar > 0 else 0.0

    BTR = Bstar * fstar / fstarTR if fstarTR > 0 else Bstar
    TTR = 2.0e5
    ReflCoef = 0.5

    for _ in range(50):
        ratioalpha = ReflCoef * (1.0 + ReflCoef) / (1.0 + ReflCoef**2) ** 1.5 * math.sqrt(2.0)
        sqbrac = (ratioalpha * Qstar * mh**2 / (rhostar**0.25 * lammax)) if lammax > 0 else 0.0
        rhoTR = sqbrac ** (4.0 / 7.0) * fstar ** (2.0 * (1.0 - theta) / 7.0)
        QTR = (Qstar * ratioalpha * (rhoTR / rhostar) ** 0.25 * math.sqrt(BTR / Bstar)) if Bstar > 0 else 0.0
        VATR = BTR / math.sqrt(4.0 * c.PI * rhoTR) if rhoTR > 0 else 1e30
        ReflCoefnew = abs((VATR - vesc) / (VATR + vesc))
        ReflCoef = math.sqrt(ReflCoefnew * ReflCoef)

    PTR = rhoTR * kb * TTR / (mu * mh)
    FluxhTR = QTR * rstar * h
    FATR = FA * fstar / fstarTR if fstarTR > 0 else FA
    if FluxhTR > FATR:
        FluxhTR = FATR

    crad = 1.4e6 * math.sqrt(lammax / lammaxsol) if lammaxsol > 0 else 0.0
    Fluxcond = crad * PTR
    if Fluxcond > 0.9:
        Fluxcond = 0.9

    mdotstarhot = ((4.0 * c.PI * rstar**2 * fstarTR / vesc2) * (FluxhTR - Fluxcond)) if vesc2 > 0 else 0.0

    uTR_hot = (
        mdotstarhot / (rhoTR * 4.0 * c.PI * rstar**2 * fstarTR) if (rhoTR > 0 and rstar > 0 and fstarTR > 0) else 0.0
    )
    MATR = uTR_hot / VATR if VATR > 0 else 0.0

    # Cold wind
    rcrit = rstar * 7.0 / (4.0 * (1.0 + (vper / vesc2) ** 2))
    ucrit = math.sqrt(c.G * M / (2.0 * rcrit)) if rcrit > 0 else 0.0
    Bcrit = Bstar * (rstar / rcrit) ** 2 * fstar if rcrit > 0 else 0.0
    vpercrit = 2.0 * ucrit

    if Bcrit > 0 and rcrit > 0 and vpercrit > 0:
        rhocrit = (
            4.0
            * c.PI
            * (rhostar * vper**2 * VA * fstar * 4.0 * c.PI * rstar**2 / (vpercrit**2 * Bcrit * 4.0 * c.PI * rcrit**2))
            ** 2
        )
        for _ in range(50):
            VAcrit = Bcrit / math.sqrt(4.0 * c.PI * rhocrit) if rhocrit > 0 else 1e30
            MAcrit = ucrit / VAcrit if VAcrit > 0 else 0.0
            macfac = (1.0 + 3.0 * MAcrit) / (1.0 + MAcrit)
            vpercrit = 2.0 * ucrit / math.sqrt(macfac) if macfac > 0 else 1.0
            rhocrit = (
                rhostar * vper**2 * VA * fstar * rstar**2 / (vpercrit**2 * VAcrit * rcrit**2 * (1.0 + MAcrit) * 2.0)
            )
    else:
        rhocrit = 0.0

    mdotstarcold = 4.0 * c.PI * rcrit**2 * ucrit * rhocrit

    mdotstar = mdotstarhot * math.exp(-4.0 * MATR**2) + mdotstarcold

    return mdotstar, Bstar * fstar, Ro, fstar
