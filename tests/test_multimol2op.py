from unittest.mock import mock_open, patch
from pathlib import Path
import shutil
from pautodock.multimol2op import (
    read_molname,
    get_mol2_name,
    split_mol2
)

def test_read_molname():
    mock_file_content = """@<TRIPOS>MOLECULE
TestMolecule
1 0 0
@<TRIPOS>ATOM
1 C 0.0 0.0 0.0 0.0 A
"""
    with patch("builtins.open", mock_open(read_data=mock_file_content)):
        molname = read_molname("test.mol2")
        assert molname == "TestMolecule"

def test_read_molname_no_molecule():
    mock_file_content = """@<TRIPOS>ATOM
1 C 0.0 0.0 0.0 0.0 A
"""
    with patch("builtins.open", mock_open(read_data=mock_file_content)):
        molname = read_molname("test.mol2")
        assert molname is None

def test_split_mol2():
    dir = Path('example_splitted')
    dir.mkdir(exist_ok=True)
    res = split_mol2('data/example.mol2', str(dir.absolute()))
    shutil.rmtree(dir.absolute())
    assert len(res) == 2
    for item in res:
        assert 'Water.mol2' in item or 'Methane.mol2' in item
    