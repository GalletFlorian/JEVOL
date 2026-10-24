"""Tests for the main evolution loop."""

import constants as c
import numpy as np
import pytest
from evolution import compute_taudec, run_evolution
from stellar_model import compute_structure, find_indices


class TestComputeTaudec:
    def test_returns_init(self):
        result = compute_taudec(1e-6, 1e-6, 1e10)
        assert result == 1e10

    def test_constant_regardless_of_velocity(self):
        r1 = compute_taudec(1e-6, 2e-6, 1e10)
        r2 = compute_taudec(5e-6, 1e-6, 1e10)
        assert r1 == r2


class TestEvolutionSynthetic:
    def test_runs_without_error(self, synthetic_config, synthetic_structure):
        n = 100
        struct = synthetic_structure
        k2rad = struct["k2rad"]
        tstar = struct["tstar"]
        tdisk = synthetic_config["tdiskparam"]

        indices = find_indices(tstar, k2rad, tdisk, 1.0, n)
        results = run_evolution(struct, synthetic_config, indices)

        assert "wconv" in results
        assert "wrad" in results

    def test_wconv_positive(self, synthetic_config, synthetic_structure):
        n = 100
        struct = synthetic_structure
        indices = find_indices(struct["tstar"], struct["k2rad"], synthetic_config["tdiskparam"], 1.0, n)
        results = run_evolution(struct, synthetic_config, indices)

        wconv = results["wconv"]
        assert np.all(wconv >= 0)

    def test_wrad_nonnegative(self, synthetic_config, synthetic_structure):
        n = 100
        struct = synthetic_structure
        indices = find_indices(struct["tstar"], struct["k2rad"], synthetic_config["tdiskparam"], 1.0, n)
        results = run_evolution(struct, synthetic_config, indices)

        wrad = results["wrad"]
        assert np.all(wrad >= 0)

    def test_amconv_decreases_over_time(self, synthetic_config, synthetic_structure):
        n = 100
        struct = synthetic_structure
        indices = find_indices(struct["tstar"], struct["k2rad"], synthetic_config["tdiskparam"], 1.0, n)
        results = run_evolution(struct, synthetic_config, indices)

        amconv = results["amconv"]
        nonzero = amconv[amconv > 0]
        if len(nonzero) > 10:
            late = nonzero[-10:]
            assert np.all(late <= nonzero[0] * 1.1)

    def test_dj_sign(self, synthetic_config, synthetic_structure):
        n = 100
        struct = synthetic_structure
        indices = find_indices(struct["tstar"], struct["k2rad"], synthetic_config["tdiskparam"], 1.0, n)
        results = run_evolution(struct, synthetic_config, indices)

        wconv = results["wconv"]
        wrad = results["wrad"]
        dj = results["dj"]

        for j in range(n):
            if wconv[j] < wrad[j]:
                assert dj[j] >= 0
            elif wconv[j] > wrad[j]:
                assert dj[j] <= 0
            else:
                assert dj[j] == 0.0


class TestEvolutionRealModel:
    @pytest.fixture(scope="class")
    def real_evolution(self):
        import os

        from parameters import build_config
        from stellar_model import load_stellar_model

        param_file = os.path.join(os.path.dirname(__file__), "..", "..", "rotevoldec.par")
        model_dir = os.path.join(os.path.dirname(__file__), "..", "..", "Model")

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

    def test_completes(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        assert len(results["wconv"]) == n

    def test_wconv_positive(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        assert np.all(results["wconv"] >= 0)

    def test_wrad_nonnegative(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        assert np.all(results["wrad"] >= 0)

    def test_final_wconv_reasonable(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        final = results["wconv"][n - 1]
        wsol_ratio = final / c.WSOL
        assert 0.001 < wsol_ratio < 1000.0

    def test_taudec_constant(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        taudec = results["taudec"]
        assert np.all(taudec == taudec[0])

    def test_total_am_decreases(self, real_evolution):
        results, struct, config, indices, n = real_evolution
        amtot = results["amconv"] + results["amrad"]
        nonzero_am = amtot[amtot > 0]
        if len(nonzero_am) > 20:
            first_nonzero = nonzero_am[0]
            last_val = nonzero_am[-1]
            assert last_val < first_nonzero * 2.0
