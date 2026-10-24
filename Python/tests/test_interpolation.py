"""Tests for sub-stepping interpolation routines."""

import constants as c
import numpy as np
from interpolation import interpol3, interpolA, interpolB


def _make_config():
    import constants as const

    n = 100
    return {
        "mstar": 1.0,
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
        "taudecinit": 561.0 * const.YEAR,
        "brklaw": 0,
        "_rstar": np.linspace(3.0, 1.0, n) * const.RSOL,
        "_lum": np.linspace(10.0, 1.0, n),
        "_teff": np.linspace(4000.0, 5800.0, n),
        "_k2conv": np.full(n, 0.2),
        "_k2rad": np.zeros(n),
        "_mrad": np.zeros(n),
        "_rrad": np.zeros(n),
        "_tstar": np.linspace(1e4, 1e9, n) * const.YEAR,
    }


class TestInterpolA:
    def test_small_change_passthrough(self):
        config = _make_config()
        r = c.RSOL
        0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        wdown = 0.01 * w0
        wup = 0.01 * w0
        result, total_down = interpolA(wdown, wup, w0, w0, r, 1.0, 5778.0, config["_tstar"], 50, 1.0, config)
        assert abs(result - (w0 + wup - wdown)) / w0 < 0.01

    def test_large_change_substeps(self):
        config = _make_config()
        r = c.RSOL
        w0 = 3.0 * c.WSOL
        wdown = 2.0 * w0
        wup = 0.5 * w0
        result, total_down = interpolA(wdown, wup, w0, w0, r, 1.0, 5778.0, config["_tstar"], 50, 1.0, config)
        assert result > 0
        assert total_down > 0


class TestInterpolB:
    def test_small_change_passthrough(self):
        config = _make_config()
        n = 100
        Iri = np.full(n, 0.01 * 1.0 * c.MSOL * c.RSOL**2)
        Ici = np.full(n, 0.2 * 1.0 * c.MSOL * c.RSOL**2)
        rri = np.full(n, 0.1 * c.RSOL)
        mri = np.full(n, 0.01)
        wr0 = 3.0 * c.WSOL
        wc0 = 3.0 * c.WSOL
        Juprad = 0.001 * Iri[50] * c.WSOL
        result, taudec = interpolB(
            Juprad, 0.0, 0.0, 0.0, wr0, wc0, Iri, Ici, rri, mri, config["_tstar"], 50, 1.0, config
        )
        assert result > 0


class TestInterpol3:
    def test_runs_without_error(self):
        config = _make_config()
        wc0 = 3.0 * c.WSOL
        wr0 = 3.0 * c.WSOL
        result_wc, result_wr = interpol3(0.0, 0.0, wc0, wr0, 0.0, 0.0, 0.0, 0.0, 0.0, 50, 1.0, config)
        assert isinstance(result_wc, float)
        assert isinstance(result_wr, float)

    def test_positive_results(self):
        config = _make_config()
        wc0 = 3.0 * c.WSOL
        wr0 = 3.0 * c.WSOL
        result_wc, result_wr = interpol3(0.0, 0.0, wc0, wr0, 0.0, 0.0, 0.0, 0.0, 0.0, 50, 1.0, config)
        assert result_wc > 0
        assert result_wr >= 0
