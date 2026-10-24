"""Tests for parameter file parsing."""

import os
import tempfile

import constants as c
from parameters import parse_parameters


class TestParseParameters:
    def test_parses_default_file(self, raw_params):
        assert "pinit" in raw_params
        assert "masse" in raw_params

    def test_d_notation_converted(self, param_file):
        params = parse_parameters(param_file)
        for _key, val in params.items():
            if isinstance(val, float):
                assert "d" not in str(val).lower() or val == 0.0

    def test_skips_comments(self):
        content = (
            "c this is a comment line 1\n"
            "c this is a comment line 2\n"
            "c line 3\n"
            "c line 4\n"
            "c line 5\n"
            "c line 6\n"
            "c line 7\n"
            "c line 8\n"
            "c line 9\n"
            "c line 10\n"
            "pinit              10.0d0\n"
            "masse              1.0d0\n"
        )
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".par", delete=False)
        f.write(content)
        f.flush()
        f.close()
        try:
            params = parse_parameters(f.name)
            assert params["pinit"] == 10.0
            assert params["masse"] == 1.0
        finally:
            os.unlink(f.name)

    def test_model_file_detected(self):
        content = (
            "c line 1\n"
            "c line 2\n"
            "c line 3\n"
            "c line 4\n"
            "c line 5\n"
            "c line 6\n"
            "c line 7\n"
            "c line 8\n"
            "c line 9\n"
            "c line 10\n"
            "modele evolutif    mod05moa16Star.dat\n"
        )
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".par", delete=False)
        f.write(content)
        f.flush()
        f.close()
        try:
            params = parse_parameters(f.name)
            assert params.get("model") == "mod05moa16Star.dat"
        finally:
            os.unlink(f.name)


class TestBuildConfig:
    def test_ksk_scaled(self, config):
        raw = parse_parameters(os.path.join(os.path.dirname(__file__), "..", "..", "rotevoldec.par"))
        raw_ksk = raw.get("ksk", 0.0)
        assert config["ksk"] == raw_ksk * c.KSK_SCALE

    def test_kmm_scaled(self, config):
        raw = parse_parameters(os.path.join(os.path.dirname(__file__), "..", "..", "rotevoldec.par"))
        raw_kmm = raw.get("kmm", 0.0)
        assert config["kmm"] == raw_kmm * c.KMM_SCALE

    def test_taudec_in_seconds(self, config):
        assert config["taudecinit"] > 1e10

    def test_tdisk_in_seconds(self, config):
        assert config["tdiskparam"] > 1e13

    def test_brklaw_is_int(self, config):
        assert isinstance(config["brklaw"], int)

    def test_model_file_string(self, config):
        assert isinstance(config["model_file"], str)
        assert ".dat" in config["model_file"]
