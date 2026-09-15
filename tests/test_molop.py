import platform
from pathlib import Path
from unittest.mock import mock_open, patch

import pytest

from pautodock.molop import (
    Molecule,
    Receptor,
    extract_coordinates,
    get_first_pose_baricentre,
    get_mol_baricentre,
    nsplit,
)


def test_default_split():
    assert nsplit("a b c") == ["a", "b", "c"]
    assert nsplit("a  b   c") == ["a", "b", "c"]


def test_custom_delimiter():
    assert nsplit("a,b,c", ",") == ["a", "b", "c"]
    assert nsplit("a,,b,,,c", ",") == ["a", "b", "c"]


def test_empty_string():
    assert nsplit("") == []


def test_all_delimiters():
    assert nsplit("   ") == []
    assert nsplit(",,,", ",") == []


def test_mixed_content():
    assert nsplit("a, b , c", ",") == ["a", " b ", " c"]


def test_leading_trailing_delimiters():
    assert nsplit(" a b c ") == ["a", "b", "c"]
    assert nsplit(",a,b,c,", ",") == ["a", "b", "c"]


def test_multicharacter_delimiter():
    assert nsplit("a||b||c", "||") == ["a", "b", "c"]
    assert nsplit("a||b||||c", "||") == ["a", "b", "c"]


def test_no_splits():
    assert nsplit("abc") == ["abc"]
    assert nsplit("abc", ",") == ["abc"]


def test_whitespace_delimiter():
    assert nsplit("a\tb\nc") == ["a", "b", "c"]
    assert nsplit("a\t\nb  c") == ["a", "b", "c"]


@pytest.mark.parametrize(
    "input_string, delimiter, expected",
    [
        ("a b c", None, ["a", "b", "c"]),
        ("a,b,c", ",", ["a", "b", "c"]),
        ("", None, []),
        ("   ", None, []),
        ("a, b , c", ",", ["a", " b ", " c"]),
        (" a b c ", None, ["a", "b", "c"]),
        ("a||b||c", "||", ["a", "b", "c"]),
        ("abc", None, ["abc"]),
        ("a\tb\nc", None, ["a", "b", "c"]),
    ],
)
def test_nsplit_parametrized(input_string, delimiter, expected):
    assert nsplit(input_string, delimiter) == expected


def test_extract_coordinates():
    line = "ATOM      1  N   ALA A   1      -0.525   1.362   0.000"
    result = extract_coordinates(line, "pdbqt")
    assert result == [-0.525, 1.362, 0.000]


def test_extract_coordinates_no_match():
    line = "REMARK This is a comment"
    result = extract_coordinates(line, "pdbqt")
    assert result is None


def test_get_mol_baricentre_pdb():
    input_ligand = "data/3EML/ligand.pdb"
    result = get_mol_baricentre(input_ligand)
    expected = [-9.06364, -7.1446, 55.8626]
    for i, val in enumerate(result):
        assert abs(val - expected[i]) < 1e-5


def test_get_mol_baricentre_pdbqt():
    input_ligand = "data/3EML/ligand.pdbqt"
    result = get_mol_baricentre(input_ligand)
    expected = [-9.06364, -7.1446, 55.8626]
    for i, val in enumerate(result):
        assert abs(val - expected[i]) < 1e-5


@pytest.fixture
def receptor():
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        return Receptor("test.pdb")


def test_receptor_init(receptor):
    assert receptor.receptor == "test.pdb"
    assert receptor.obabel_path == "/usr/bin/"


def test_receptor_topdbqt(tmp_path):
    rec = tmp_path / "test.pdb"
    rec.write_text(
        "ATOM      1  CA  ALA A   1      "
        + "%8.3f%8.3f%8.3f"
        % (
            1.0,
            2.0,
            3.0,
        )
        + "  1.00  0.00           C\n"
    )
    pdbqt = tmp_path / "test.pdbqt"
    pdbqt.write_text(
        "REMARK  Name = test\n"
        "ROOT\n"
        "ATOM      1  CA  ALA A   1      "
        + "%8.3f%8.3f%8.3f" % (1.0, 2.0, 3.0)
        + "  1.00  0.00           C\n"
        "ENDROOT\n"
        "BRANCH   1   2\n"
        "HETATM    2  O   HOH A 200      "
        + "%8.3f%8.3f%8.3f" % (4.0, 5.0, 6.0)
        + "  1.00  0.00    -0.834 OA\n"
        "ENDBRANCH\n"
        "TORSDOF 0\n"
    )
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        receptor = Receptor(str(rec))
        with patch("subprocess.run") as mock_run:
            result = receptor.topdbqt()
            mock_run.assert_called_once()
    assert result.endswith("test.pdbqt")
    # only the ATOM record must survive: the water molecule, the
    # REMARK lines and the ligand-style tree records are stripped,
    # like in the prepare_receptor4.py output
    lines = pdbqt.read_text().splitlines()
    assert lines == [
        "ATOM      1  CA  ALA A   1      "
        + "%8.3f%8.3f%8.3f" % (1.0, 2.0, 3.0)
        + "  1.00  0.00           C",
    ]


def test_receptor_mgltools_method(tmp_path):
    rec = tmp_path / "test.pdb"
    rec.write_text("")
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        receptor = Receptor(str(rec), method="mgltools")
        assert receptor.method == "mgltools"
        with patch("pathlib.Path.exists", return_value=True):
            with patch("subprocess.call") as mock_call:
                result = receptor.topdbqt()
                mock_call.assert_called_once()
                assert result.endswith("test.pdbqt")


def test_receptor_unknown_method():
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        with pytest.raises(ValueError, match="Unknown receptor preparation method"):
            Receptor("test.pdb", method="foo")


@pytest.fixture
def molecule():
    system = platform.system()
    patch_value = None
    if system == "Linux":
        patch_value = "/usr/bin/"
    elif system == "Darwin":
        patch_value = "/opt/homebrew/bin/"

    with patch("pautodock.fileutils.get_bin_path", return_value=patch_value):
        return Molecule("test.mol2")


def test_molecule_init(molecule):
    assert molecule.molecule == "test.mol2"
    system = platform.system()
    patch_value = None
    if system == "Linux":
        patch_value = "/usr/bin/"
    elif system == "Darwin":
        patch_value = "/opt/homebrew/bin/"
    assert molecule.obabel_path == patch_value


def test_molecule_topdbqt(molecule):
    with patch("subprocess.call") as mock_call:
        with patch("builtins.open", mock_open()):
            result = molecule.topdbqt()
            mock_call.assert_called_once()
            assert result.endswith("test.pdbqt")


def test_molecule_topdbqt_with_translation(tmp_path):
    mol2 = tmp_path / "test.mol2"
    mol2.write_text("@<TRIPOS>MOLECULE\nTestMol\n")
    prefix = "ATOM      1  N10 ZMA A 401    "
    pdbqt = tmp_path / "test.pdbqt"
    pdbqt.write_text(
        prefix
        + "%8.3f%8.3f%8.3f" % (-9.420, -9.544, 56.644)
        + "  0.00  0.00    +0.000 NA\n"
        + prefix
        + "%8.3f%8.3f%8.3f" % (-8.953, -8.593, 55.842)
        + "  0.00  0.00    +0.000 C\n"
        + "REMARK unchanged line\n"
    )
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        with patch("subprocess.call") as mock_call:
            mol = Molecule(str(mol2))
            result = mol.topdbqt([1.0, 2.0, 3.0])
            mock_call.assert_called_once()
    # the molecule baricentre must land on the target center
    cc = get_mol_baricentre(result)
    assert abs(cc[0] - 1.0) < 1e-3
    assert abs(cc[1] - 2.0) < 1e-3
    assert abs(cc[2] - 3.0) < 1e-3
    assert "REMARK unchanged line" in Path(result).read_text()


def test_extract_coordinates_spec_aligned():
    # 30-char prefix, coordinates in the spec fields 31-54 (1-based)
    line = (
        "ATOM      1  C1  LIG A   1    "
        + "%8.3f%8.3f%8.3f" % (-123.456, 7.123, -0.5)
        + "  0.00  0.00    +0.000 C"
    )
    assert extract_coordinates(line, "pdbqt") == [-123.456, 7.123, -0.5]
    assert extract_coordinates(line, "pdb") == [-123.456, 7.123, -0.5]


def test_get_mol_baricentre_pdb_atom_records(tmp_path):
    pdb = tmp_path / "prot.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1    "
        + "%8.3f%8.3f%8.3f" % (1.0, 2.0, 3.0)
        + "  1.00  0.00           C\n"
        "ATOM      2  CA  ALA A   2    "
        + "%8.3f%8.3f%8.3f" % (3.0, 4.0, 5.0)
        + "  1.00  0.00           C\n"
    )
    assert get_mol_baricentre(str(pdb)) == [2.0, 3.0, 4.0]


def test_get_first_pose_baricentre(tmp_path):
    pdbqt = tmp_path / "poses.pdbqt"
    prefix = "ATOM      1  N10 ZMA A 401    "
    pdbqt.write_text(
        "MODEL 1\n"
        + prefix
        + "%8.3f%8.3f%8.3f" % (1.0, 1.0, 1.0)
        + "  0.00  0.00    +0.000 NA\n"
        + "ENDMDL\n"
        "MODEL 2\n"
        + prefix
        + "%8.3f%8.3f%8.3f" % (50.0, 50.0, 50.0)
        + "  0.00  0.00    +0.000 NA\n"
        + "ENDMDL\n"
    )
    assert get_first_pose_baricentre(str(pdbqt)) == [1.0, 1.0, 1.0]


def test_molecule_topdbqt_with_ph(tmp_path):
    mol2 = tmp_path / "test.mol2"
    mol2.write_text("@<TRIPOS>MOLECULE\nTestMol\n")
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        with patch("subprocess.call") as mock_call:
            mol = Molecule(str(mol2))
            mol.topdbqt(ph=7.4)
            cmd = mock_call.call_args[0][0][0]
            assert "-p 7.4" in cmd


def test_molecule_topdbqt_without_ph(tmp_path):
    mol2 = tmp_path / "test.mol2"
    mol2.write_text("@<TRIPOS>MOLECULE\nTestMol\n")
    with patch("pautodock.molop.get_bin_path", return_value="/usr/bin/"):
        with patch("subprocess.call") as mock_call:
            mol = Molecule(str(mol2))
            mol.topdbqt()
            cmd = mock_call.call_args[0][0][0]
            assert " -p " not in cmd
