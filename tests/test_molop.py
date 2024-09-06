import pytest
import platform
from unittest.mock import patch, mock_open
from pautodock.molop import(
    nsplit,
    extract_coordinates,
    get_mol_baricentre,
    Receptor,
    Molecule
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

@pytest.mark.parametrize("input_string, delimiter, expected", [
    ("a b c", None, ["a", "b", "c"]),
    ("a,b,c", ",", ["a", "b", "c"]),
    ("", None, []),
    ("   ", None, []),
    ("a, b , c", ",", ["a", " b ", " c"]),
    (" a b c ", None, ["a", "b", "c"]),
    ("a||b||c", "||", ["a", "b", "c"]),
    ("abc", None, ["abc"]),
    ("a\tb\nc", None, ["a", "b", "c"]),
])
def test_nsplit_parametrized(input_string, delimiter, expected):
    assert nsplit(input_string, delimiter) == expected
                             
def test_extract_coordinates():
    line = "ATOM      1  N   ALA A   1      -0.525   1.362   0.000"
    result = extract_coordinates(line, 'pdbqt')
    assert result == [-0.525, 1.362, 0.000]

def test_extract_coordinates_no_match():
    line = "REMARK This is a comment"
    result = extract_coordinates(line, 'pdbqt')
    assert result is None

def test_get_mol_baricentre_pdb():
    input_ligand = 'data/3EML/ligand.pdb'
    result = get_mol_baricentre(input_ligand)
    expected = [-9.06364, -7.1446, 55.8626]
    for i, val in enumerate(result):
        assert abs(val-expected[i]) < 1e-5

def test_get_mol_baricentre_pdbqt():
    input_ligand = 'data/3EML/ligand.pdbqt'
    result = get_mol_baricentre(input_ligand)
    expected = [-9.06364, -7.1446, 55.8626]
    for i, val in enumerate(result):
        assert abs(val-expected[i]) < 1e-5

@pytest.fixture
def receptor():
    return Receptor("test.pdb", "/path/to/mgl")

def test_receptor_init(receptor):
    assert receptor.receptor == "test.pdb"
    assert receptor.mglpath == "/path/to/mgl"

def test_receptor_topdbqt(receptor):
    with patch('subprocess.call') as mock_call:
        result = receptor.topdbqt()
        mock_call.assert_called_once()
        assert result.endswith("test.pdbqt")

@pytest.fixture
def molecule():
    system = platform.system()
    patch_value = None
    if system == 'Linux': 
        patch_value='/usr/bin/'
    elif system == 'Darwin':
        patch_value='/opt/homebrew/bin/'

    with patch('pautodock.fileutils.get_bin_path', return_value=patch_value):
        return Molecule("test.mol2", "/path/to/mgl")

def test_molecule_init(molecule):
    assert molecule.molecule == "test.mol2"
    assert molecule.mglpath == "/path/to/mgl"
    system = platform.system()
    patch_value = None
    if system == 'Linux': 
        patch_value='/usr/bin/'
    elif system == 'Darwin':
        patch_value='/opt/homebrew/bin/'
    assert molecule.obabel_path == patch_value

def test_molecule_topdbqt(molecule):
    with patch('subprocess.call') as mock_call:
        with patch('builtins.open', mock_open()):
            result = molecule.topdbqt()
            mock_call.assert_called_once()
            assert result.endswith("test.pdbqt")

def test_molecule_topdbqt_with_translation(molecule):
    mock_file_content = """
ATOM      1  N10 ZMA A 401      -9.420  -9.544  56.644  0.00  0.00    +0.000 NA
ATOM      2  C11 ZMA A 401      -8.953  -8.593  55.842  0.00  0.00    +0.000 C 
ATOM      3  N12 ZMA A 401      -8.585  -8.857  54.578  0.00  0.00    +0.000 NA
    """
    with patch('subprocess.call'):
        with patch('builtins.open', mock_open(read_data=mock_file_content)) as mock_file:
            result = molecule.topdbqt([1.0, 1.0, 1.0])
            assert result.endswith("test.pdbqt")
            mock_file.assert_called
