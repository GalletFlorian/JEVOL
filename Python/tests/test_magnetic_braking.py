"""Tests for magnetic braking calculations."""

import constants as c
import numpy as np
from magnetic_braking import magnetic_braking


def _make_config(
    brklaw=0, ksk=0.0, kmm=0.0, ksc=0.0, kmp=0.1e1, K1MP=8.5, K2MP=0.0506, m=0.22, K3=2.8, K4=0.2, mvic=0.11, K=1.0
):
    import constants as const

    return {
        "mstar": 1.0,
        "ksk": ksk,
        "kmm": kmm,
        "ksc": ksc,
        "kmp": kmp,
        "K": K,
        "K1MP": K1MP,
        "K2MP": K2MP,
        "m": m,
        "K3": K3,
        "K4": K4,
        "mvic": mvic,
        "taudecinit": 561.0 * const.YEAR,
        "brklaw": brklaw,
        "_rstar": np.ones(100) * const.RSOL,
        "_lum": np.ones(100),
        "_teff": np.ones(100) * 5778.0,
        "_k2conv": np.ones(100) * 0.2,
        "_k2rad": np.zeros(100),
        "_mrad": np.zeros(100),
        "_rrad": np.zeros(100),
        "_tstar": np.linspace(1e4, 1e9, 100) * const.YEAR,
    }


class TestSkumanich:
    def test_braking_positive(self):
        config = _make_config(ksk=1.5e-12 * c.KSK_SCALE)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0

    def test_braking_zero_dt(self):
        config = _make_config(ksk=1.5e-12 * c.KSK_SCALE)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        result = magnetic_braking(w0, r, 1.0, 0.0, 0.0, ic, config, 5778.0)
        assert result == 0.0

    def test_braking_scales_with_w3(self):
        config = _make_config(ksk=1.5e-12 * c.KSK_SCALE)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        dt = 1e6 * c.YEAR
        w0_a = 2.0 * c.WSOL
        w0_b = 4.0 * c.WSOL
        res_a = magnetic_braking(w0_a, r, 1.0, dt, 0.0, ic, config, 5778.0)
        res_b = magnetic_braking(w0_b, r, 1.0, dt, 0.0, ic, config, 5778.0)
        ratio = res_b / res_a
        assert abs(ratio - 8.0) < 0.1

    def test_braking_scales_with_1_ic(self):
        config = _make_config(ksk=1.5e-12 * c.KSK_SCALE)
        r = c.RSOL
        ic1 = 0.2 * 1.0 * c.MSOL * r**2
        ic2 = 0.4 * 1.0 * c.MSOL * r**2
        dt = 1e6 * c.YEAR
        w0 = 3.0 * c.WSOL
        res1 = magnetic_braking(w0, r, 1.0, dt, 0.0, ic1, config, 5778.0)
        res2 = magnetic_braking(w0, r, 1.0, dt, 0.0, ic2, config, 5778.0)
        assert abs(res1 / res2 - 2.0) < 0.01


class TestMayorMermilliod:
    def test_braking_positive(self):
        config = _make_config(kmm=1.25e-11 * c.KMM_SCALE)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0

    def test_braking_scales_with_w2(self):
        config = _make_config(kmm=1.25e-11 * c.KMM_SCALE)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        dt = 1e6 * c.YEAR
        w0_a = 2.0 * c.WSOL
        w0_b = 4.0 * c.WSOL
        res_a = magnetic_braking(w0_a, r, 1.0, dt, 0.0, ic, config, 5778.0)
        res_b = magnetic_braking(w0_b, r, 1.0, dt, 0.0, ic, config, 5778.0)
        ratio = res_b / res_a
        assert abs(ratio - 4.0) < 0.1


class TestMatt2012:
    def test_braking_positive(self):
        config = _make_config(brklaw=0, kmp=0.1e1)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0

    def test_braking_zero_dt(self):
        config = _make_config(brklaw=0, kmp=0.1e1)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        result = magnetic_braking(w0, r, 1.0, 0.0, 0.0, ic, config, 5778.0)
        assert result == 0.0


class TestReville2015:
    def test_braking_positive(self):
        config = _make_config(brklaw=1, kmp=0.1e1)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0


class TestMatt2015:
    def test_braking_positive(self):
        config = _make_config(brklaw=3, kmp=0.1e1)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0


class TestReinersMohanty:
    def test_braking_positive(self):
        config = _make_config(ksk=0.0, kmm=0.0, ksc=0.0, kmp=0.0)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        w0 = 3.0 * c.WSOL
        dt = 1e6 * c.YEAR
        result = magnetic_braking(w0, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert result > 0

    def test_sub_vs_super_critical(self):
        config = _make_config(ksk=0.0, kmm=0.0, ksc=0.0, kmp=0.0)
        r = c.RSOL
        ic = 0.2 * 1.0 * c.MSOL * r**2
        dt = 1e6 * c.YEAR
        Wcr = 8.56e-6
        w_sub = 0.5 * Wcr
        w_sup = 2.0 * Wcr
        res_sub = magnetic_braking(w_sub, r, 1.0, dt, 0.0, ic, config, 5778.0)
        res_sup = magnetic_braking(w_sup, r, 1.0, dt, 0.0, ic, config, 5778.0)
        assert res_sub > 0
        assert res_sup > 0
        assert res_sup > res_sub
