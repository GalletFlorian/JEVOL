# JEVOL - Stellar Rotational Evolution Simulator

JEVOL simulates the rotational evolution of low-mass stars (0.5-1.5 Msol) from
birth to main sequence, including:
- Core-envelope decoupling
- Magnetized wind braking
- Star-disk coupling during pre-main-sequence evolution
- Hot and cold stellar wind mass loss

It is a Python translation of the Fortran code `rotevoldec.f` by F. Gallet.

## Quick start

```bash
# From the project root (JEVOL/)
cd Python
python main.py

# Or from anywhere
python Python/main.py
```

Output is written to `Python_opencode/output/`. See [USER_GUIDE.md](USER_GUIDE.md)
for parameter configuration and output file descriptions.

## Repository structure

```
JEVOL/
  rotevoldec.par              # Simulation parameter file
  rotevoldec.f                # Original Fortran code (reference)
  magnbrak.f                  # Fortran magnetic braking
  mdot.f                      # Fortran mass loss
  interlin.f                  # Fortran interpolation
  interspline.f               # Fortran spline interpolation
  Model/                      # Stellar evolution model files
    mod10moa16Star.dat        # 1.0 Msol model
    mod05moa16Star.dat        # 0.5 Msol model
  Python_opencode/            # Python implementation
    main.py                   # Entry point
    parameters.py             # Parameter file parser
    stellar_model.py          # Model loader and structure computation
    evolution.py              # Main evolution loop (3 phases)
    magnetic_braking.py       # 7 braking laws
    mass_loss.py              # Cranmer hot+cold wind model
    interpolation.py          # Sub-stepping routines
    output_writer.py          # Output file generation
    constants.py              # Physical constants (CGS)
    conftest.py               # Pytest fixtures
    tests/                    # 81 unit tests
    output/                   # Simulation output
  doc/                        # This documentation
    README.md                 # Overview (this file)
    USER_GUIDE.md             # How to run and configure simulations
    DEVELOPER_GUIDE.md        # Code architecture and Fortran mapping
    PHYSICS.md                # Physical models and equations
```

## Performance

The Python code runs in ~3-5 seconds for a 1.0 Msol model (1334 timesteps).
This is achieved through:
- Reduced sub-stepping in Phase III interpolation (ninter=20 vs Fortran's 2000)
- Eliminated redundant dict copies in hot loops
- Class-scoped test fixtures (test suite runs in ~8s)

## Requirements

- Python 3.10+
- numpy

## Testing

```bash
cd Python
python -m pytest tests/ -v
```

81 tests covering constants, parameter parsing, stellar model loading,
evolution, interpolation, magnetic braking, mass loss, and integration.

## Known differences from Fortran

The Python code is a faithful translation with a few intentional improvements:
- `1/3` integer division bug fixed (Fortran had `(1/3)=0` in Reiners & Mohanty)
- `r**(3/2)` integer division fixed (Fortran had `r**1= r`)
- `Fluxcond` correctly uses thermal pressure PTR instead of temperature TTR
- K1MP reduction for massive cores implemented (Fortran line 554-556)

See [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md) for full comparison details.
