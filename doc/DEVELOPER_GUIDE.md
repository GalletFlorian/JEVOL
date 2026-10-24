# Developer Guide

## Code architecture

The simulation follows a simple pipeline:

```
main.py
  -> parameters.build_config()        # Parse rotevoldec.par
  -> stellar_model.load_stellar_model()  # Read model file
  -> stellar_model.compute_structure()   # Compute Ienv, Icore, etc.
  -> stellar_model.find_indices()        # Find phase boundaries
  -> evolution.run_evolution()           # Main loop
  -> output_writer.write_outputs()       # Write .dat files
```

### Module responsibilities

| Module | Lines | Role |
|--------|-------|------|
| `main.py` | 67 | Entry point, orchestrates pipeline |
| `parameters.py` | 77 | Parses `rotevoldec.par`, builds config dict |
| `stellar_model.py` | 167 | Loads model files, computes structure, finds indices |
| `evolution.py` | 248 | Main loop: Phase I, II.A/II.B, III |
| `magnetic_braking.py` | 135 | 7 braking laws (dispatch by coefficients) |
| `mass_loss.py` | 156 | Cranmer hot+cold wind model, Bstar, Rossby number |
| `interpolation.py` | 227 | Sub-stepping: interpolA, interpolB, interpol3 |
| `output_writer.py` | 164 | Writes all .dat output files |
| `constants.py` | 30 | Physical constants in CGS |

## Evolution phases

The evolution is divided into 3 main phases based on two critical times:
- **tdisk**: when the protoplanetary disk disappears
- **tdec**: when a radiative core appears (k2rad crosses 1e-9)

### Phase I: Disk-locked (steps 1 to min(tdisk, tdec))
Envelope rotation is fixed at `winit = 2*pi / (pinit * day)`.
Core has not formed yet (wrad = 0).

### Phase II: Between disk loss and core formation
Two sub-cases depending on the order of tdisk and tdec:

**Case A** (tdisk < tdec, typical for low-mass):
Disk disappears first. Envelope spins down under magnetic braking
while core is still absent. Uses `interpolA` sub-stepping.

**Case B** (tdec < tdisk, typical for higher mass):
Core forms first. Core spins up while envelope remains disk-locked.
Uses `interpolB` sub-stepping. At the first step, wrad is forced
equal to wconv.

### Phase III: Both disk gone and core formed (steps max(tdisk,tdec)+1 to end)
Both core and envelope evolve. Envelope loses angular momentum via
magnetic braking; core-envelope coupling transfers angular momentum.
Uses `interpol3` sub-stepping.

## Fortran-to-Python mapping

Each Python module maps to a Fortran source:

| Python | Fortran | Notes |
|--------|---------|-------|
| `evolution.py` | `rotevoldec.f` lines 438-850 | Main loop |
| `magnetic_braking.py` | `magnbrak.f` | Braking laws |
| `mass_loss.py` | `mdot.f` | Mass loss + magnetic field |
| `interpolation.py` | `interlin.f` | Sub-stepping interpolation |
| `parameters.py` | `rotevoldec.f` lines 100-170 | Parameter parsing |
| `stellar_model.py` | `rotevoldec.f` lines 170-400 | Structure + indices |
| `constants.py` | `rotevoldec.f` lines 150-170 | Physical constants |
| `output_writer.py` | `rotevoldec.f` lines 859-1005 | Output generation |

### Key Fortran line references

| Fortran lines | Python location | Description |
|---------------|-----------------|-------------|
| magnbrak.f:130-133 | magnetic_braking.py:78-84 | Matt 2012 braking formula |
| magnbrak.f:151-178 | magnetic_braking.py:100-121 | Matt 2015 braking formula |
| mdot.f:213 | mass_loss.py:113 | PTR = rhoTR * kb * TTR / (mu*mh) |
| mdot.f:226 | mass_loss.py:120 | Fluxcond = crad * PTR |
| rotevoldec.f:554-556 | evolution.py:113-115 | K1MP reduction for massive cores |
| rotevoldec.f:634-682 | evolution.py:145-163 | Phase II.B deltaJ computation |
| interlin.f:92 | interpolation.py:87-93 | interpolA total braking sum |
| interlin.f:246 | interpolation.py:242-243 | interpol3 taudec forcing |

## Interpolation sub-stepping

The three interpolation routines divide large timesteps into smaller sub-steps
when velocity changes exceed a threshold:

- **interpolA**: Phase II.A (envelope evolving, no core). Dynamic ninter
  based on `Wtest / (0.1 * |W0|)`, clamped to [10, 200].
- **interpolB**: Phase II.B (core evolving, disk-locked envelope). Same
  dynamic ninter formula.
- **interpol3**: Phase III (both evolving). Fixed ninter=20 (Fortran uses
  2000, but compiled code handles this).

Each sub-step:
1. Interpolates stellar structure (Ienv, Icore, rrad, mrad) linearly
2. Interpolates luminosity and Teff
3. Calls `magnetic_braking` for the sub-step
4. Updates Wconv and Wrad with structure changes + braking + coupling

## Configuration dict

The `config` dict is the central data carrier. Keys include:

```python
{
    # Stellar parameters
    "mstar": float,         # Stellar mass (Msol)
    "pinit": float,         # Initial period (days)
    "model_file": str,      # Model filename

    # Braking parameters
    "ksk": float,           # Skumanich coefficient (scaled)
    "kmm": float,           # Mayor-Mermilliod coefficient (scaled)
    "ksc": float,           # Schatzman coefficient (scaled)
    "kmp": float,           # Matt efficiency
    "brklaw": int,          # Braking law selector
    "K": float,             # Matt 2012 parameter
    "K1MP": float,          # Matt 2012 magnetic field scaling
    "K2MP": float,          # Matt 2012 escape velocity scaling
    "m": float,             # Magnetic field exponent
    "K3": float,            # Reville constant
    "K4": float,            # Reville critical velocity
    "mvic": float,          # Reville velocity exponent

    # Timescales
    "taudecinit": float,    # Coupling timescale (seconds)
    "tdiskparam": float,    # Disk lifetime (seconds)

    # Stellar evolution arrays (set during run_evolution)
    "_lum": ndarray,        # Luminosity array
    "_teff": ndarray,       # Effective temperature array
    "_rstar": ndarray,      # Radius array (cm)
    "_tstar": ndarray,      # Time array (seconds)
    "_k2conv": ndarray,     # Convective k2
    "_k2rad": ndarray,      # Radiative k2
    "_mrad": ndarray,       # Core mass
    "_rrad": ndarray,       # Core radius (cm)
}
```

## Testing

81 tests organized by module:

| File | Tests | Coverage |
|------|-------|----------|
| `test_constants.py` | 13 | All physical constants |
| `test_parameters.py` | 8 | Parsing, config building |
| `test_stellar_model.py` | 11 | Loading, structure, indices |
| `test_evolution.py` | 11 | Taudec, synthetic + real model |
| `test_interpolation.py` | 5 | interpolA, interpolB, interpol3 |
| `test_magnetic_braking.py` | 11 | All 7 braking laws |
| `test_mass_loss.py` | 7 | Mdot, Bstar, Rossby, filling factor |
| `test_integration.py` | 8 | Full simulation, output files |

Fixtures use `scope="class"` to avoid re-running the full simulation
for each test method.

```bash
python -m pytest tests/ -v
```

## Performance considerations

The main bottleneck is `mdot_calculation` (~80% of runtime), called
~25K times (1262 Phase III steps x 20 sub-steps). Each call runs two
50-iteration convergence loops (reflection coefficient, cold wind density).

Key optimizations applied:
1. Removed `dict(config)` copies in interpolation hot loops (252K eliminated)
2. Reduced `ninter` from 200 to 20 in `interpol3` (max error 0.46%)
3. Class-scoped pytest fixtures

To tune accuracy vs speed, adjust `ninter` in `interpolation.py:171`.
Values of 10-50 are reasonable for the default 1 Msol model.

## Adding a new braking law

1. Add the formula to `magnetic_braking.py` as a new `elif` branch in
   `magnetic_braking()`.
2. Add any new parameters to `parameters.py:build_config()`.
3. Add the parameter to `rotevoldec.par`.
4. Add unit tests in `tests/test_magnetic_braking.py`.
5. Update the case numbering in the `if/elif` chain (currently cases 1-5).

## Adding a new stellar model

1. Place the model file in `Model/`.
2. Format: 1-line header, then columns `idx Lum Teff tstar rstar k2conv k2rad Mrad Rrad`.
3. Set `modele evolutif <filename>` in `rotevoldec.par`.
4. Adjust `masse` to match the model mass.
