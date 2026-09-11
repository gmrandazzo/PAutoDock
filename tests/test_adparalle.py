import platform
import tempfile
from unittest.mock import patch

import pytest

from pautodock.adparallel import ADParallel


@pytest.fixture
def ad_parallel():
    system = platform.system()
    patch_value = None
    if system == "Linux":
        patch_value = "/usr/bin/"
    elif system == "Darwin":
        patch_value = "/opt/homebrew/bin/"

    with patch("pautodock.adparallel.get_bin_path", return_value=patch_value):
        with patch("pautodock.adparallel.install_mgltools"):
            receptor = "path/to/receptor.pdb"
            ligand = "path/to/ligand.mol2"
            db = "path/to/database.mol2"
            wpath = tempfile.mkdtemp()
            return ADParallel(receptor, ligand, db, wpath)


def test_init(ad_parallel):
    assert ad_parallel.receptor == "path/to/receptor.pdb"
    assert ad_parallel.ligand == "path/to/ligand.mol2"
    assert ad_parallel.db == "path/to/database.mol2"
    assert isinstance(ad_parallel.wpath, str)
    assert ad_parallel.gsize_x == 30
    assert ad_parallel.gsize_y == 30
    assert ad_parallel.gsize_z == 30
    assert ad_parallel.speed == "fast"
    assert ad_parallel.atd is True
    assert ad_parallel.vina is True


def test_read_atom_types(ad_parallel, tmp_path):
    test_pdb = tmp_path / "test.pdb"
    test_pdb.write_text(
        "ATOM   3511  CA  TYR A1161     -10.160  10.285   7.878  1.00 56.11     0.191 C\n"  # noqa: E501
        "ATOM   3512  C   TYR A1161     -10.552  10.249   9.353  1.00 64.26     0.275 C\n"  # noqa: E501
        "ATOM   3513  O   TYR A1161      -9.767   9.828  10.203  1.00 63.26    -0.268 OA\n"  # noqa: E501
    )
    atypes_str, atlst = ad_parallel.read_atom_types(str(test_pdb))
    assert atypes_str == "OA C" or atypes_str == "C OA"
    assert atlst == ["OA", "C"] or atlst == ["C", "OA"]


def test_write_vina_param_files(ad_parallel, tmp_path):
    cc = [1.0, 2.0, 3.0]
    ss = [30, 30, 30]

    vina_conf_path = ad_parallel.write_vina_param_files(str(tmp_path), cc, ss)

    assert vina_conf_path.exists()
    assert vina_conf_path.name == "vina_conf.txt"

    with open(vina_conf_path, "r") as f:
        content = f.read()
        assert "center_x = 1.0000" in content
        assert "center_y = 2.0000" in content
        assert "center_z = 3.0000" in content
        assert "size_x = 30" in content
        assert "size_y = 30" in content
        assert "size_z = 30" in content


def test_write_autodock_param_files_speed_modes(ad_parallel, tmp_path):
    rec_pdbqt = tmp_path / "rec.pdbqt"
    rec_pdbqt.write_text(
        "ATOM      1  CA  TYR A1161     -10.160  10.285   7.878  1.00 56.11     0.191 C\n"  # noqa: E501
    )
    mol_pdbqt = tmp_path / "mol.pdbqt"
    mol_pdbqt.write_text(
        "ATOM      1  CA  TYR A1161     -10.160  10.285   7.878  1.00 56.11     0.191 C\n"  # noqa: E501
    )
    for speed, evals in [
        ("fast", "ga_num_evals 250000 "),
        ("normal", "ga_num_evals 2500000 "),
        ("thorough", "ga_num_evals 25000000 "),
    ]:
        ad_parallel.speed = speed
        _, dpf_path = ad_parallel.write_autodock_param_files(
            str(tmp_path), str(rec_pdbqt), mol_pdbqt.name, [0.0, 0.0, 0.0]
        )
        with open(dpf_path, "r") as f:
            assert evals in f.read()


def test_read_vina_output(ad_parallel, tmp_path):
    vina_out = tmp_path / "vina_out.txt"
    vina_out.write_text(
        "mode |   affinity | dist from best mode\n"
        "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
        "-----+------------+----------+----------\n"
        "   1       -9.317          0          0\n"
        "   2       -8.823      3.083      9.394\n"
        "   3       -8.819      2.787      9.676\n"
    )

    avg, min_val, max_val = ad_parallel.read_vina_output(str(vina_out))

    assert abs(-8.986 - avg) < 1e-4
    assert (-9.317 - min_val) < 1e-4
    assert (-8.819 - max_val) < 1e-4
