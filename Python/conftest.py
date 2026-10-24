"""Shared test fixtures."""

import os
import sys

import pytest

# Add Python to path for sibling imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import constants as c
from parameters import build_config, parse_parameters
from stellar_model import compute_structure, load_stellar_model

BASE_DIR = os.path.join(os.path.dirname(__file__), "..")
PARAM_FILE = os.path.join(BASE_DIR, "rotevoldec.par")
MODEL_DIR = os.path.join(BASE_DIR, "Model")
DEFAULT_MODEL = os.path.join(MODEL_DIR, "mod05moa16Star.dat")


@pytest.fixture
def param_file():
    return PARAM_FILE


@pytest.fixture
def model_dir():
    return MODEL_DIR


@pytest.fixture
def raw_params():
    if os.path.exists(PARAM_FILE):
        return parse_parameters(PARAM_FILE)
    pytest.skip("rotevoldec.par not found")


@pytest.fixture
def config():
    if os.path.exists(PARAM_FILE):
        return build_config(PARAM_FILE)
    pytest.skip("rotevoldec.par not found")


@pytest.fixture
def loaded_model():
    if os.path.exists(DEFAULT_MODEL):
        return load_stellar_model(DEFAULT_MODEL)
    pytest.skip("Default model file not found")


@pytest.fixture
def structure(loaded_model, config):
    n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
    struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, config["mstar"])
    config["_lum"] = lum
    config["_teff"] = teff
    return struct, config, n


@pytest.fixture
def synthetic_config():
    """Minimal config for unit tests without file I/O."""
    return {
        "mstar": 1.0,
        "pinit": 10.0,
        "ksk": 0.0,
        "kmm": 0.0,
        "ksc": 0.0,
        "kmp": 0.1e1,
        "K": 1.0,
        "K1MP": 8.5,
        "K2MP": 0.0506,
        "m": 0.22,
        "K3": 2.8,
        "K4": 0.2,
        "mvic": 0.11,
        "taudecinit": 561.0 * c.YEAR,
        "tdiskparam": 4.0e6 * c.YEAR,
        "brklaw": 0,
    }


@pytest.fixture
def synthetic_structure():
    """Synthetic stellar structure for unit tests."""
    import numpy as np

    n = 100
    tstar_yr = np.linspace(1e4, 1e9, n)
    rstar_rsol = np.linspace(3.0, 1.0, n)
    k2conv_arr = np.full(n, 0.2)
    k2rad_arr = np.where(tstar_yr > 5e7, 0.01, 0.0)
    mrad_arr = np.where(tstar_yr > 5e7, 0.01, 0.0)
    rrad_rel = np.where(tstar_yr > 5e7, 0.1, 0.0)
    lum = np.linspace(10.0, 1.0, n)
    teff = np.linspace(4000.0, 5800.0, n)

    struct = compute_structure(tstar_yr, rstar_rsol, k2conv_arr, k2rad_arr, mrad_arr, rrad_rel, 1.0)
    struct["lum"] = lum
    struct["teff"] = teff
    return struct
