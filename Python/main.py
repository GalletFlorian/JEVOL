"""
JEVOL - Stellar Rotational Evolution Simulator (Python version)
Entry point.

Can be run as:
    python Python/main.py        (from project root)
    python -m Python.main        (from project root)
"""

import os
import sys

# Add script directory to path for sibling imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from evolution import run_evolution
from output_writer import write_outputs
from parameters import build_config
from stellar_model import compute_structure, find_indices, load_stellar_model


def main(param_file=None, model_dir=None, output_dir=None):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.normpath(os.path.join(script_dir, ".."))

    if param_file is None:
        param_file = os.path.join(project_root, "rotevoldec.par")
    if model_dir is None:
        model_dir = os.path.join(project_root, "Model")
    if output_dir is None:
        output_dir = os.path.join(script_dir, "output")

    config = build_config(param_file)
    print(f"Stellar mass: {config['mstar']} Msol")
    print(f"Initial period: {config['pinit']} days")
    print(f"Braking law: {config['brklaw']}")

    model_file = os.path.join(model_dir, config["model_file"])
    n, tstar_yr, rstar_rsol, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(model_file)
    print(f"Loaded {n} timesteps from {config['model_file']}")

    struct = compute_structure(tstar_yr, rstar_rsol, k2conv, k2rad, mrad, rrad, config["mstar"])

    config["_lum"] = lum
    config["_teff"] = teff

    tdisk_sec = config["tdiskparam"]
    indices = find_indices(struct["tstar"], k2rad, tdisk_sec, config["mstar"], n)

    print(f"Disk disappears at step: {indices[0]}")
    print(f"Core decoupling at step: {indices[1]}")
    print(f"Main sequence at step: {indices[2]}")

    results = run_evolution(struct, config, indices)

    write_outputs(results, struct, config, indices, output_dir)
    print(f"Output written to {output_dir}/")

    print("Done.")
    return results, struct, config, indices


if __name__ == "__main__":
    main()
