"""
Output writer for JEVOL simulation results.
"""

import os

import constants as c
import numpy as np


def write_outputs(results, struct, config, indices, output_dir="output"):
    """
    Write key output .dat files.

    Parameters:
        results: dict from evolution.run_evolution
        struct: dict from stellar_model.compute_structure
        config: dict from parameters.build_config
        indices: tuple from stellar_model.find_indices
        output_dir: directory to write output files
    """
    os.makedirs(output_dir, exist_ok=True)

    n = len(struct["tstar"])
    tstar = struct["tstar"]
    rstar = struct["rstar"]
    Ienv = struct["Ienv"]
    Icore = struct["Icore"]
    wbreakup = struct["wbreakup"]
    struct["vbreakup"]

    wconv = results["wconv"]
    wrad = results["wrad"]
    dj = results["dj"]
    amconv = results["amconv"]
    amrad = results["amrad"]
    taudec = results["taudec"]

    config["mstar"]
    ndebut = indices[4]

    time_myr = tstar / (c.YEAR * 1.0e6)

    pconv = np.where(wconv > 0, c.DPI / np.where(wconv > 0, wconv, 1.0) / c.DAY, 0.0)
    prad = np.where(wrad > 0, c.DPI / np.where(wrad > 0, wrad, 1.0) / c.DAY, 0.0)

    vconv = rstar * wconv / 1.0e5
    vrad = rstar * wrad / 1.0e5

    wconv_wsol = wconv / c.WSOL
    wrad_wsol = wrad / c.WSOL

    jtot = (amconv + amrad) / c.AMSOL

    files = {}

    def _open(name):
        path = os.path.join(output_dir, name)
        f = open(path, "w")
        files[name] = f
        return f

    def _write(name, fmt, *cols):
        if name not in files:
            _open(name)
        rows = zip(*[col[ndebut:n] for col in cols], strict=False)
        for row in rows:
            files[name].write(fmt.format(*row) + "\n")

    fmt_time_val = "{:14.8f} {:12.6e}\n"

    # Write each output file
    with _open("wconvevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], wconv_wsol[i]))

    with _open("wradevol.dat") as f:
        for i in range(ndebut, n):
            if wrad[i] != 0.0:
                f.write(fmt_time_val.format(time_myr[i], wrad_wsol[i]))

    with _open("pconvevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], pconv[i]))

    with _open("pradevol.dat") as f:
        for i in range(ndebut, n):
            if wrad[i] != 0.0:
                f.write(fmt_time_val.format(time_myr[i], prad[i]))

    with _open("vconvevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], vconv[i]))

    with _open("vradevol.dat") as f:
        for i in range(ndebut, n):
            if wrad[i] != 0.0:
                f.write(fmt_time_val.format(time_myr[i], vrad[i]))

    with _open("jevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], jtot[i]))

    with _open("jconvevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], amconv[i] / c.AMSOL))

    with _open("jradevol.dat") as f:
        for i in range(ndebut, n):
            f.write(fmt_time_val.format(time_myr[i], amrad[i] / c.AMSOL))

    with _open("djevol.dat") as f:
        for i in range(ndebut, n):
            if taudec[i] > 0:
                f.write(fmt_time_val.format(time_myr[i], dj[i] / taudec[i]))

    with _open("wbreakup.dat") as f:
        for i in range(ndebut, n):
            f.write(f"{time_myr[i]:10.4f} {wbreakup[i]:12.6e}\n")

    with _open("wconvwbu.dat") as f:
        for i in range(ndebut, n):
            if wbreakup[i] > 0:
                f.write(fmt_time_val.format(time_myr[i], wconv[i] / wbreakup[i]))

    with _open("wradwbu.dat") as f:
        for i in range(ndebut, n):
            if wbreakup[i] > 0:
                f.write(fmt_time_val.format(time_myr[i], wrad[i] / wbreakup[i]))

    with _open("taudec.dat") as f:
        for i in range(ndebut, n):
            f.write(f"{time_myr[i]:10.4f} {taudec[i] / (c.YEAR * 1e6):12.6e}\n")

    with _open("Iconvdtevol.dat") as f:
        for i in range(ndebut, n):
            f.write(f"{time_myr[i]:10.4f} {Ienv[i] / c.ISOL:12.6e}\n")

    with _open("Iraddtevol.dat") as f:
        for i in range(ndebut, n):
            f.write(f"{time_myr[i]:10.4f} {Icore[i] / c.ISOL:12.6e}\n")

    _write_parameter_summary(config, output_dir)

    for f in files.values():
        f.close()


def _write_parameter_summary(config, output_dir):
    path = os.path.join(output_dir, "ReadmeParam.dat")
    with open(path, "w") as f:
        f.write("JEVOL Python - Parameter Summary\n")
        f.write("=" * 40 + "\n")
        f.write(f"Pinit = {config['pinit']} days\n")
        f.write(f"Mstar = {config['mstar']} Msol\n")
        f.write(f"Model = {config['model_file']}\n")
        f.write(f"Braking law = {config['brklaw']}\n")
        f.write(f"K1MP = {config['K1MP']}\n")
        f.write(f"K2MP = {config['K2MP']}\n")
        f.write(f"m = {config['m']}\n")
        f.write(f"Taudec = {config['taudecinit'] / c.YEAR} yr\n")
        f.write(f"Tdisk = {config['tdiskparam'] / c.YEAR} yr\n")
