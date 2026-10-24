# Physical Models

## Overview

JEVOL computes the rotational evolution of a low-mass star by solving
the angular momentum equation for two coupled zones (convective envelope
and radiative core):

```
dJ_env/dt = -T_wind - T_coupling
dJ_core/dt = +T_coupling
```

where:
- `T_wind` is the magnetic braking torque (wind removes angular momentum)
- `T_coupling` transfers angular momentum between core and envelope
  proportional to their velocity difference and inversely proportional
  to the coupling timescale taudec

## Evolution phases

### Phase I: Disk-coupled
The star is locked to the protoplanetary disk at the initial period `pinit`.
The envelope rotation is constant: `wconv = winit`.
No radiative core exists yet.

Duration: from birth until min(tdisk, tdec).

### Phase II.A: Disk gone, no core (typical for 1 Msol)
The disk has disappeared but the radiative core has not yet formed.
The envelope spins down under magnetic braking alone.
Sub-stepping via `interpolA`.

### Phase II.B: Core formed, disk still present
The radiative core has formed but the disk is still present.
The envelope remains disk-locked (`wconv = winit`).
The core spins up as it gains mass from the contracting interior.
At the first step, `wrad` is forced equal to `wconv`.
Sub-stepping via `interpolB`.

### Phase III: Fully decoupled
Both core and envelope evolve. The envelope loses angular momentum
through magnetic braking. The core-envelope coupling transfers
angular momentum proportionally to their velocity difference.

Sub-stepping via `interpol3` with ninter=20 sub-steps per timestep.

## Magnetic braking laws

The braking law determines how much angular momentum the stellar wind
carries away. Five families are implemented:

### Case 1: Skumanich (1972)

```
dW/dt = -ksk * W^3 * (R/Rsol)^0.5 * M^(-0.5) / I
```

Classic braking law. Angular velocity decays as `W ~ t^(-1/2)`.
Active when `ksk != 0`.

### Case 2: Mayor & Mermilliod (1991)

```
dV/dt = -kmm * V^2 * (R/Rsol)^0.5 * M^(-0.5) / I
```

Surface velocity braking. A saturation regime (`ksat`) applies
when `W > wsat = 30 * Wsol`, switching to linear braking.
Active when `kmm != 0`.

### Case 3: Schatzman (1962)

```
dW/dt = -ksc * W^(7/3) * (R/Rsol)^0.5 * M^(-0.5) / I
```

Steeper than Skumanich. Active when `ksc != 0`.

### Case 4: Matt et al. (2012) / Reville et al. (2015) / Matt et al. (2015)

Active when `kmp != 0` and `ksk = kmm = 0`. Requires mass loss rate
computation (see Mass loss section).

**brklaw=0 (Matt 2012):**
```
T_wind = K1MP^2 * kmp * Bstar^(4m) * R^(5m+2) * W * Mdot^(1-2m) * dt
         / (I * (K2MP^2 * 2*G*M + K*W^2*R^3)^m)
```

**brklaw=1,2 (Reville 2015):**
```
T_wind = K3^2 * Bstar^(4*mvic) * R^(5*mvic+2) * W * Mdot^(1-2*mvic)
         * (4pi)^(4*mvic) * dt
         / (I * (2*G*M + (2/K4^2)*W^2*R^3)^mvic)
```

**brklaw=3 (Matt 2015):**
Uses convective turnover time `tcz` to determine the braking regime.
Two branches based on Rossby number Ro vs Rosun/chi:
- Fast rotators (Ro < Rosun/chi): linear in W
- Slow rotators (Ro > Rosun/chi): power law `W^(p+1)` with p=2.3

### Case 5: Reiners & Mohanty (2012)

```
if W < Wcr:
    dW/dt = C * (W/Wcr)^4 * W * (R^16/M^2)^(1/3) * dt / I
else:
    dW/dt = C * W * (R^16/M^2)^(1/3) * dt / I
```

Used when all k coefficients are zero. Critical velocity Wcr = 8.56e-6 rad/s.

## Mass loss model (Cranmer et al. 2007)

The mass loss rate determines the wind torque in Case 4 braking laws.

### Hot wind

Based on energy flux deposition in the transition region:

```
Mdot_hot = (4*pi*R^2*fstarTR/Vesc^2) * min(Flux_hTR, FATR - Fluxcond)
```

where:
- `Flux_hTR` = coronal heating flux at the transition region
- `Fluxcond` = conductive flux (uses thermal pressure PTR)
- `FATR` = maximum available flux

### Cold wind

Based on Alfen-wave-driven wind at the critical point:

```
Mdot_cold = 4*pi*rcrit^2 * ucrit * rhocrit
```

where rcrit is the critical point radius, ucrit the sound speed there,
and rhocrit the density determined iteratively (50 iterations).

### Combined mass loss

```
Mdot = Mdot_hot * exp(-4*MATR^2) + Mdot_cold
```

where MATR = uTR/VATR is the Alfen Mach number at the transition region.

### Magnetic field

```
Bstar = 1.13 * sqrt(8*pi*rhostar*kb*Teff/(mu*mh))
```

The factor 1.13 accounts for flux tube expansion from the photosphere
to the coronal base.

### Filling factor

```
f = 0.55 / (1 + (Ro/Rosol / 0.16)^2.3)^1.22
```

where Ro = Omega/tauc is the Rossby number and tauc is the convective
turnover time.

## Core-envelope coupling

The angular momentum exchange between core and envelope:

```
dJ_coupling = (Icore * Ienv / (Icore + Ienv)) * (Wrad - Wconv) * dt / taudec
```

This is applied symmetrically:
- Core loses `dJ_coupling` (slows down if Wrad > Wconv)
- Envelope gains `dJ_coupling` (speeds up)

The coupling timescale `taudec` is constant (forced to the initial value)
despite a Spada et al. formula being implemented in the code.

## Sub-stepping interpolation

When the velocity change in a single timestep exceeds ~10% of the current
velocity, the timestep is subdivided into `ninter` sub-steps. During each
sub-step, the stellar structure (Ienv, Icore, Rrad, Mrad) is linearly
interpolated.

Three interpolation routines handle different evolution phases:
- `interpolA`: envelope only (Phase II.A)
- `interpolB`: core only (Phase II.B)
- `interpol3`: both core and envelope (Phase III)

The number of sub-steps is:
- Dynamic in interpolA/interpolB: `ninter = clamp(Wtest/(0.1*|W0|), 10, 200)`
- Fixed in interpol3: `ninter = 20`

## Units

All internal computations use CGS units:
- Length: cm
- Mass: g
- Time: s
- Angular velocity: rad/s
- Magnetic field: Gauss
- Energy: erg

Solar reference values:
- WSOL = 2.87e-6 rad/s
- RSOL = 6.9599e10 cm
- MSOL = 1.989e33 g
- LSOL = 3.826e33 erg/s
- AMSOL = 1.63e48 g*cm^2/s
- ISOL = 6.411e53 g*cm^2
