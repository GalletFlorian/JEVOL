"""Tests for stellar model loading and structure computation."""

import os

import constants as c
import numpy as np
import pytest
from stellar_model import compute_structure, find_indices, load_stellar_model

MODEL_FILE = os.path.join(os.path.dirname(__file__), "..", "..", "Model", "mod05moa16Star.dat")


class TestLoadStellarModel:
    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_loads_nonzero_steps(self):
        n, *_ = load_stellar_model(MODEL_FILE)
        assert n > 0

    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_tstar_positive(self):
        n, tstar, rstar, *_ = load_stellar_model(MODEL_FILE)
        assert np.all(tstar[:n] > 0)

    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_rstar_positive(self):
        n, tstar, rstar, *_ = load_stellar_model(MODEL_FILE)
        assert np.all(rstar[:n] > 0)

    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_k2conv_nonnegative(self):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(MODEL_FILE)
        assert np.all(k2conv >= 0)

    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_k2rad_nonnegative(self):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(MODEL_FILE)
        assert np.all(k2rad >= 0)

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_stellar_model("nonexistent_file.dat")


class TestComputeStructure:
    def test_tstar_in_seconds(self, loaded_model):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 1.0)
        assert np.all(struct["tstar"] > 1e10)

    def test_rstar_in_cm(self, loaded_model):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 1.0)
        assert np.all(struct["rstar"] > 1e10)

    def test_iconv_nonnegative(self, loaded_model):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 1.0)
        assert np.all(struct["Ienv"] >= 0)

    def test_wbreakup_positive(self, loaded_model):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 1.0)
        assert np.all(struct["wbreakup"] > 0)

    def test_vbreakup_positive(self, loaded_model):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = loaded_model
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 1.0)
        assert np.all(struct["vbreakup"] > 0)


class TestFindIndices:
    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_index_disk_in_range(self):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(MODEL_FILE)
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 0.5)
        tdisk = 4.0e6 * c.YEAR
        tstar_sec = struct["tstar"]
        idx = find_indices(tstar_sec, k2rad, tdisk, 0.5, n)
        assert 0 <= idx[0] < n

    @pytest.mark.skipif(not os.path.exists(MODEL_FILE), reason="Model file not found")
    def test_index_ms_in_range(self):
        n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(MODEL_FILE)
        struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, 0.5)
        tdisk = 4.0e6 * c.YEAR
        tstar_sec = struct["tstar"]
        idx = find_indices(tstar_sec, k2rad, tdisk, 0.5, n)
        assert 0 <= idx[2] < n
