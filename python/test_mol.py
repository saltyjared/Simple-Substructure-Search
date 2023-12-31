"""Tests for Molecule object"""

from molecule import Molecule
import pytest


@pytest.mark.parametrize(
    ("filepath"),
    [
        pytest.param("sdf/tnt.sdf", id="tnt"),
        pytest.param("sdf/vanillin.sdf", id="vanillin"),
    ],
)
def test_Molecule_instantiation(filepath):
    """
    Tests that Molecule instances are properly constructed from filepaths to
    .sdf files.
    """
    mol = Molecule(filepath)
    assert isinstance(mol, Molecule)
    assert hasattr(mol, "element_dict")
    assert hasattr(mol, "bond_tuples")
    assert hasattr(mol, "graph")


@pytest.mark.parametrize(
    ("left, right, expected"),
    [
        pytest.param("sdf/tnt.sdf", "sdf/tnt.sdf", True, id="is_equal"),
        pytest.param("sdf/tnt.sdf", "sdf/vanillin.sdf", False, id="is_not_equal"),
    ],
)
def test_Molecule_equivalence(left, right, expected):
    """
    Tests the equivalence operator on different Molecule instances based on
    their molecular fingerprint.
    """
    mol1 = Molecule(left)
    mol2 = Molecule(right)
    assert (mol1 == mol2) == expected


@pytest.mark.parametrize(
    ("molecule, substruct, expected"),
    [
        pytest.param("sdf/tnt.sdf", "sdf/tnt.sdf", True, id="has_itself"),
        pytest.param("sdf/tnt.sdf", "sdf/no2.sdf", True, id="has_no2"),
        pytest.param("sdf/tnt.sdf", "sdf/ethylene.sdf", True, id="has_C=C"),
        pytest.param("sdf/tnt.sdf", "sdf/acetylene.sdf", False, id="no_C#C"),
        pytest.param("sdf/tnt.sdf", "sdf/vanillin.sdf", False, id="no_vanillin"),
    ],
)
def test_check_substructure_sdf(molecule, substruct, expected):
    """
    Tests whether a molecule has a given substructure using its .sdf file.
    """
    mol = Molecule(molecule)
    assert mol.check_substructure_sdf(Molecule(substruct)) == expected


@pytest.mark.parametrize(
    ("molecule, substruct, expected"),
    [
        pytest.param("sdf/tnt.sdf", "ON=O", True, id="has_no2"),
        pytest.param("sdf/tnt.sdf", "C=C", True, id="has_C=C"),
        pytest.param("sdf/tnt.sdf", "C#C", False, id="no_C#C"),
    ],
)
def test_check_substructure_smiles(molecule, substruct, expected):
    """
    Tests whether a molecule has a given substructure using its SMILES path.
    """
    mol = Molecule(molecule)
    assert mol.check_substructure_smiles(substruct) == expected


def test_equivalence_error():
    """
    Tests that when using the equivalence operator, both objects must be a
    Molecule.
    """
    mol1 = Molecule("sdf/tnt.sdf")
    mol2 = "Hello, world!"

    with pytest.raises(TypeError):
        mol1 == mol2


def test_check_sub_sdf_error():
    """
    Tests that when using the check_substructure_sdf method, the input must be
    a Molecule.
    """
    mol1 = Molecule("sdf/tnt.sdf")
    mol2 = "Hello, world!"

    with pytest.raises(TypeError):
        mol1.check_substructure_sdf(mol2)


def test_check_sub_smiles_error():
    """
    Tests that when using the check_substructure_sdf method, the input must be
    a string.
    """
    mol1 = Molecule("sdf/tnt.sdf")
    mol2 = Molecule("sdf/tnt.sdf")

    with pytest.raises(TypeError):
        mol1.check_substructure_smiles(mol2)
