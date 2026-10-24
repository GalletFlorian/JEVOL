"""
Parameter file parser for rotevoldec.par.
"""


def parse_parameters(param_file):
    """
    Parse rotevoldec.par parameter file.

    Skips the first 10 comment/header lines, then reads key-value pairs
    where keys are descriptive labels and values are the actual parameters.
    """
    params = {}

    with open(param_file) as f:
        lines = f.readlines()

    for line in lines[10:]:
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        if parts[0].lower() == "parameter":
            key = parts[1].lower()
            raw_val = parts[2].replace("d", "e").replace("D", "e") if len(parts) > 2 else parts[1]
        else:
            key = parts[0].lower()
            raw_val = parts[1].replace("d", "e").replace("D", "e")

        try:
            params[key] = float(raw_val)
        except ValueError:
            for part in parts[1:]:
                if ".dat" in part.lower():
                    params["model"] = part
                    break

    return params


def build_config(param_file):
    """
    Parse parameter file and return a config dict with all simulation parameters,
    with constant multiplications and unit conversions already applied.
    """
    import constants as c

    raw = parse_parameters(param_file)

    config = {}
    config["pinit"] = raw.get("pinit", 11.0)
    config["mstar"] = raw.get("masse", 1.0)
    config["model_file"] = raw.get("model", "mod05moa16Star.dat")

    config["ksk"] = raw.get("ksk", 0.0) * c.KSK_SCALE
    config["kmm"] = raw.get("kmm", 0.0) * c.KMM_SCALE
    config["ksc"] = raw.get("ksc", 0.0) * c.KSC_SCALE
    config["kmp"] = raw.get("kmp", 0.0)
    config["K"] = raw.get("k", 1.0)
    config["K1MP"] = raw.get("k1mp", 1.0)
    config["K2MP"] = raw.get("k2mp", 1.0)
    config["m"] = raw.get("m", 0.22)
    config["K3"] = raw.get("k3", 2.8)
    config["K4"] = raw.get("k4", 0.2)
    config["mvic"] = raw.get("mvic", 0.11)

    config["taudecinit"] = raw.get("taudec", 561.0) * c.YEAR
    config["tdiskparam"] = raw.get("tdisk", 4.0e6) * c.YEAR
    config["brklaw"] = int(raw.get("brklaw", 0))

    return config
