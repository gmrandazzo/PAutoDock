import sys
from unittest.mock import patch

import pytest

from pautodock.__main__ import parse_arguments


def _argv(*extra):
    return [
        "pautodock",
        "--receptor",
        "r.pdb",
        "--wdir",
        "w",
        "--db",
        "d.mol2",
        "--cx",
        "0",
        "--cy",
        "0",
        "--cz",
        "0",
    ] + list(extra)


def test_ph_default_is_none():
    with patch.object(sys, "argv", _argv()):
        assert parse_arguments().ph is None


def test_ph_option():
    with patch.object(sys, "argv", _argv("--ph", "7.4")):
        assert parse_arguments().ph == 7.4


def test_ph_out_of_range():
    with patch.object(sys, "argv", _argv("--ph", "20")):
        with pytest.raises(SystemExit):
            parse_arguments()
