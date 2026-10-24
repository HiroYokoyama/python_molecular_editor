"""Validate chemical project geometry before document replacement."""

import base64
import copy

import pytest
from rdkit import Chem

from moleditpy.core.project_validation import (
    validate_project_data,
    validate_raw_project_data,
)


def _structure():
    return {
        "atoms": [
            {"id": 0, "symbol": "C", "x": 0, "y": 0},
            {"id": 1, "symbol": "N", "x": 50, "y": 0},
        ],
        "bonds": [{"atom1": 0, "atom2": 1, "order": 1}],
        "next_atom_id": 2,
    }


@pytest.mark.parametrize(
    "change",
    [
        lambda s: s["atoms"].append(copy.deepcopy(s["atoms"][0])),
        lambda s: s["atoms"][0].update(id=True),
        lambda s: s["atoms"][0].update(charge=2**31),
        lambda s: s["bonds"].append(copy.deepcopy(s["bonds"][0])),
        lambda s: s["bonds"][0].update(atom2=0),
        lambda s: s["bonds"][0].update(order=4),
        lambda s: s["bonds"][0].update(order=True),
        lambda s: s["bonds"][0].update(stereo=5),
        lambda s: s.update(next_atom_id=1),
    ],
)
def test_rejects_inconsistent_2d_geometry(change):
    structure = _structure()
    change(structure)
    with pytest.raises(ValueError):
        validate_project_data({"2d_structure": structure})


def _mol():
    mol = Chem.MolFromSmiles("CN")
    conformer = Chem.Conformer(2)
    conformer.SetAtomPosition(1, (1.5, 0, 0))
    mol.AddConformer(conformer)
    return mol


def test_valid_project_and_raw_geometry():
    mol = _mol()
    validate_project_data({"2d_structure": _structure(), "3d_structure": None})
    validate_project_data(
        {
            "3d_structure": {
                "mol_binary_base64": base64.b64encode(mol.ToBinary()).decode()
            }
        }
    )
    validate_raw_project_data(
        {
            "atoms": {
                0: {"symbol": "C", "pos": (0, 0)},
                1: {"symbol": "N", "pos": (50, 0)},
            },
            "bonds": {(0, 1): {"order": 1}},
            "mol_3d": mol.ToBinary(),
        }
    )


@pytest.mark.parametrize(
    "atoms",
    [
        [{"index": 0}],
        [{"index": 0}, {"index": 0}],
        [{"index": 0}, {"index": 2}],
        [{"index": 0, "x": float("inf")}, {"index": 1}],
    ],
)
def test_rejects_inconsistent_3d_metadata(atoms):
    encoded = base64.b64encode(_mol().ToBinary()).decode()
    with pytest.raises(ValueError):
        validate_project_data(
            {"3d_structure": {"mol_binary_base64": encoded, "atoms": atoms}}
        )


def test_rejects_nonfinite_binary_conformer():
    mol = _mol()
    mol.GetConformer().SetAtomPosition(0, (float("nan"), 0, 0))
    with pytest.raises(ValueError, match="finite"):
        validate_raw_project_data({"mol_3d": mol.ToBinary()})
