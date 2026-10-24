"""Integration tests - full simulation run."""

import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import constants as c
from evolution import run_evolution
from output_writer import write_outputs
from parameters import build_config
from stellar_model import compute_structure, find_indices, load_stellar_model

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..")


@pytest.fixture(scope="class")
def full_simulation():
    """Run full simulation with default parameters."""
    param_file = os.path.join(BASE_DIR, "rotevoldec.par")
    model_dir = os.path.join(BASE_DIR, "Model")

    if not os.path.exists(param_file):
        pytest.skip("rotevoldec.par not found")

    config = build_config(param_file)
    model_file = os.path.join(model_dir, config["model_file"])

    if not os.path.exists(model_file):
        pytest.skip("Model file not found")

    n, tstar, rstar, lum, teff, k2conv, k2rad, mrad, rrad, k2 = load_stellar_model(model_file)
    struct = compute_structure(tstar, rstar, k2conv, k2rad, mrad, rrad, config["mstar"])
    config["_lum"] = lum
    config["_teff"] = teff

    tdisk = config["tdiskparam"]
    indices = find_indices(struct["tstar"], k2rad, tdisk, config["mstar"], n)
    results = run_evolution(struct, config, indices)

    return results, struct, config, indices, n


class TestFullSimulation:
    def test_completes(self, full_simulation):
        results, struct, config, indices, n = full_simulation
        assert n > 0
        assert len(results["wconv"]) == n

    def test_no_nans(self, full_simulation):
        results, _, _, _, _ = full_simulation
        assert not np.any(np.isnan(results["wconv"]))
        assert not np.any(np.isnan(results["wrad"]))
        assert not np.any(np.isnan(results["dj"]))

    def test_no_infs(self, full_simulation):
        results, _, _, _, _ = full_simulation
        assert not np.any(np.isinf(results["wconv"]))
        assert not np.any(np.isinf(results["wrad"]))

    def test_positive_velocities(self, full_simulation):
        results, _, _, _, _ = full_simulation
        assert np.all(results["wconv"] >= 0)
        assert np.all(results["wrad"] >= 0)

    def test_output_files_created(self, full_simulation):
        results, struct, config, indices, n = full_simulation
        tmpdir = tempfile.mkdtemp()
        try:
            write_outputs(results, struct, config, indices, tmpdir)
            expected_files = [
                "wconvevol.dat",
                "wradevol.dat",
                "pconvevol.dat",
                "pradevol.dat",
                "vconvevol.dat",
                "vradevol.dat",
                "jevol.dat",
                "jconvevol.dat",
                "jradevol.dat",
                "djevol.dat",
                "wbreakup.dat",
                "taudec.dat",
                "wconvwbu.dat",
                "wradwbu.dat",
                "Iconvdtevol.dat",
                "Iraddtevol.dat",
                "ReadmeParam.dat",
            ]
            for fname in expected_files:
                path = os.path.join(tmpdir, fname)
                assert os.path.exists(path), f"Missing output: {fname}"
                assert os.path.getsize(path) > 0, f"Empty output: {fname}"
        finally:
            shutil.rmtree(tmpdir)

    def test_output_content_physical(self, full_simulation):
        results, struct, config, indices, n = full_simulation
        tmpdir = tempfile.mkdtemp()
        try:
            write_outputs(results, struct, config, indices, tmpdir)

            wconv_data = np.loadtxt(os.path.join(tmpdir, "wconvevol.dat"))
            assert wconv_data.shape[1] >= 2
            assert np.all(wconv_data[:, 0] > 0)
            assert np.all(wconv_data[:, 1] > 0)

            vconv_data = np.loadtxt(os.path.join(tmpdir, "vconvevol.dat"))
            assert np.all(vconv_data[:, 1] > 0)
            assert np.all(vconv_data[:, 1] < 1000.0)

            taudec_data = np.loadtxt(os.path.join(tmpdir, "taudec.dat"))
            assert np.all(taudec_data[:, 1] > 0)

        finally:
            shutil.rmtree(tmpdir)

    def test_period_evolution(self, full_simulation):
        results, struct, config, indices, n = full_simulation
        wconv = results["wconv"]
        valid = wconv[wconv > 0]
        if len(valid) > 10:
            periods = c.DPI / valid / c.DAY
            assert np.all(periods > 0)
            assert np.all(periods < 10000)

    def test_angular_momentum_budget(self, full_simulation):
        results, struct, config, indices, n = full_simulation
        amtot = results["amconv"] + results["amrad"]
        for j in range(1, n):
            if amtot[j] > 0 and amtot[j - 1] > 0:
                assert amtot[j] <= amtot[j - 1] * 1.01
