"""Tests for mass loss rate calculation."""

import constants as c
from mass_loss import mdot_calculation


def _make_config():
    import constants as const

    return {
        "mstar": 1.0,
        "taudecinit": 561.0 * const.YEAR,
    }


class TestMdotCalculation:
    def test_solar_like_star(self):
        config = _make_config()
        r = c.RSOL
        w = c.WSOL
        lum = 1.0
        teff = 5778.0
        mdot, bstar, ro, fstar = mdot_calculation(w, r, lum, config, teff)
        assert mdot > 0
        assert bstar > 0
        assert ro > 0
        assert 0 < fstar < 1

    def test_mdot_reasonable_range(self):
        config = _make_config()
        r = c.RSOL
        w = c.WSOL
        lum = 1.0
        teff = 5778.0
        mdot, _, _, _ = mdot_calculation(w, r, lum, config, teff)
        mdot_msol_yr = mdot / c.MDOTYR
        assert 1e-16 < mdot_msol_yr < 1e-8

    def test_bstar_solar_range(self):
        config = _make_config()
        w = c.WSOL
        r = c.RSOL
        lum = 1.0
        teff = 5778.0
        _, bstar, _, _ = mdot_calculation(w, r, lum, config, teff)
        assert 0.1 < bstar < 10000.0

    def test_rossby_slow_rotator(self):
        config = _make_config()
        r = c.RSOL
        w = 0.5 * c.WSOL
        lum = 1.0
        teff = 5778.0
        _, _, ro, _ = mdot_calculation(w, r, lum, config, teff)
        assert ro > 1.0

    def test_rossby_fast_rotator(self):
        config = _make_config()
        r = c.RSOL
        w = 100.0 * c.WSOL
        lum = 1.0
        teff = 5778.0
        _, _, ro, _ = mdot_calculation(w, r, lum, config, teff)
        assert ro < 1.0

    def test_filling_factor_bounded(self):
        config = _make_config()
        r = c.RSOL
        lum = 1.0
        teff = 5778.0
        for w in [0.1 * c.WSOL, c.WSOL, 10.0 * c.WSOL, 100.0 * c.WSOL]:
            _, _, _, fstar = mdot_calculation(w, r, lum, config, teff)
            assert 0 < fstar <= 1.0, f"fstar={fstar} for w={w}"

    def test_mdot_increases_with_activity(self):
        config = _make_config()
        r = c.RSOL
        lum = 1.0
        teff = 5778.0
        mdot_slow, _, _, _ = mdot_calculation(0.5 * c.WSOL, r, lum, config, teff)
        mdot_fast, _, _, _ = mdot_calculation(50.0 * c.WSOL, r, lum, config, teff)
        assert mdot_fast > mdot_slow
