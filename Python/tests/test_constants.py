"""Tests for physical constants."""

import constants as c


def test_solar_mass():
    assert abs(c.MSOL - 1.989e33) < 1e30


def test_solar_radius():
    assert abs(c.RSOL - 6.9599e10) < 1e8


def test_solar_angular_velocity():
    assert abs(c.WSOL - 2.87e-6) < 1e-8


def test_solar_luminosity():
    assert abs(c.LSOL - 3.826e33) < 1e30


def test_solar_angular_momentum():
    assert abs(c.AMSOL - 1.63e48) < 1e45


def test_solar_moment_of_inertia():
    assert abs(c.ISOL - 6.411e53) < 1e50


def test_gravitational_constant():
    assert abs(c.G - 6.673e-8) < 1e-10


def test_day_in_seconds():
    assert abs(c.DAY - 86400.0) < 1.0


def test_year_in_seconds():
    expected = 365.0 * 86400.0
    assert abs(c.YEAR - expected) < 100.0


def test_dpi():
    import math

    assert abs(c.DPI - 2.0 * math.pi) < 1e-4


def test_pi():
    assert abs(c.PI - 3.14159) < 1e-5


def test_age_sun():
    assert c.AGESUN == 4.5e9


def test_constants_are_positive():
    assert c.MSOL > 0
    assert c.RSOL > 0
    assert c.WSOL > 0
    assert c.G > 0
    assert c.DAY > 0
    assert c.YEAR > 0
