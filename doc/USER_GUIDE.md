# User Guide

## Running a simulation

```bash
python Python/main.py
```

The simulation reads `rotevoldec.par` (parameter file) and a stellar model file
from `Model/`. Output is written to `Python_opencode/output/`.

### Custom parameters

Edit `rotevoldec.par` in the project root. The first 10 lines are comments.
After that, each line is a key-value pair:

```
pinit              8.0d0                day
masse              1.0d0                 solar mass
modele evolutif    mod10moa16Star.dat
```

Fortran `d` notation (e.g. `1.7d0`) is automatically converted to `1.7e0`.

### Parameter reference

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `pinit` | 11.0 | days | Initial rotation period (constant until disk decoupling) |
| `masse` | 1.0 | Msol | Stellar mass |
| `modele evolutif` | — | — | Stellar model filename (in `Model/`) |
| `taudec` | 561.0 | Myr | Core-envelope coupling timescale |
| `tdisk` | 4.0 | Myr | Disk lifetime |
| `brklaw` | 0 | — | Braking law selector (see below) |

#### Braking law parameters

| Parameter | Used by | Description |
|-----------|---------|-------------|
| `ksk` | Skumanich | Wind braking coefficient (w^3 law) |
| `kmm` | Mayor-Mermilliod | Wind braking coefficient (w^2 law) |
| `ksc` | Schatzman | Wind braking coefficient (w^(7/3) law) |
| `kmp` | Matt 2012+ | Wind braking efficiency |
| `parameter K` | Matt 2012 | Rotation parameter |
| `parameter K1MP` | Matt 2012 | Magnetic field scaling |
| `parameter K2MP` | Matt 2012 | Escape velocity scaling |
| `parameter m` | Matt 2012, 2015 | Magnetic field exponent |
| `parameter K3` | Reville 2015 | Braking constant |
| `parameter K4` | Reville 2015 | Critical velocity scaling |
| `parameter mvic` | Reville 2015 | Velocity exponent |

Note: Parameters prefixed with `parameter` (e.g. `parameter K1MP 1.7d0`) are
read correctly. The `parameter` keyword is stripped and the second token is
the key.

### Braking law selection (`brklaw`)

| Value | Law | Physics |
|-------|-----|---------|
| 0 | Matt et al. 2012 | B ~ rho^(1/2), dJ/dt from wind torque |
| 1 | Reville et al. 2015 | Generalized wind torque |
| 2 | Reville et al. 2015 | Same as 1 (multi-epoch weighted) |
| 3 | Matt et al. 2015 | Convective turnover time dependent |
| — | Skumanich | dW/dt ~ -W^3 (use `ksk != 0`) |
| — | Mayor-Mermilliod | dV/dt ~ -V^2 (use `kmm != 0`) |
| — | Schatzman | dW/dt ~ -W^(7/3) (use `ksc != 0`) |
| — | Reiners & Mohanty | Used when all k=0 |

The braking law is selected by which coefficient is nonzero.
`ksk`, `kmm`, `ksc` are pre-scaled by constants in the code.

### Stellar model files

Located in `Model/`. Format (after 1-line header):

```
idx  Lum  Teff  tstar(yr)  rstar(Rsol)  k2conv  k2rad  Mrad  Rrad
```

Available models:
- `mod10moa16Star.dat` — 1.0 Msol, 1334 timesteps
- `mod05moa16Star.dat` — 0.5 Msol, 400 timesteps

## Output files

All output is written to `Python_opencode/output/` as two-column ASCII files:
`time(Myr)  value`.

### Velocity evolution

| File | Quantity | Units |
|------|----------|-------|
| `wconvevol.dat` | Envelope angular velocity | Wsol |
| `wradevol.dat` | Core angular velocity | Wsol |
| `vconvevol.dat` | Envelope surface velocity | km/s |
| `vradevol.dat` | Core surface velocity | km/s |
| `wconvwbu.dat` | Envelope W / breakup W | dimensionless |
| `wradwbu.dat` | Core W / breakup W | dimensionless |
| `wbreakup.dat` | Breakup angular velocity | Wsol |

### Period evolution

| File | Quantity | Units |
|------|----------|-------|
| `pconvevol.dat` | Envelope rotation period | days |
| `pradevol.dat` | Core rotation period | days |

### Angular momentum

| File | Quantity | Units |
|------|----------|-------|
| `jevol.dat` | Total angular momentum | AMSOL |
| `jconvevol.dat` | Envelope angular momentum | AMSOL |
| `jradevol.dat` | Core angular momentum | AMSOL |
| `djevol.dat` | Core-envelope exchange rate | AMSOL/s |

### Structure

| File | Quantity | Units |
|------|----------|-------|
| `Iconvdtevol.dat` | Envelope moment of inertia | ISOL |
| `Iraddtevol.dat` | Core moment of inertia | ISOL |
| `taudec.dat` | Coupling timescale | Myr |

### Summary

| File | Description |
|------|-------------|
| `ReadmeParam.dat` | Run parameters summary |

### Plotting

A plotting script is available:

```bash
python Python/plot.py
```

Generates `Python_opencode/output/evolution_plot.png`.
