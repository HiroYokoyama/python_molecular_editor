"""Unit tests for core.stereo_check: 2D-drawn vs 3D chirality comparison."""

import ast
import os
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

import moleditpy.core as core_pkg
from moleditpy.core.molecular_data import MolecularData
from moleditpy.core.stereo_check import (
    actual_chirality,
    drawn_chirality,
    find_chirality_mismatches,
)
from moleditpy.utils.constants import ANGSTROM_PER_PIXEL

# Tetrodotoxin, PubChem CID 11174599 (fetched 2026-09-23). The IUPAC name lists
# nine stereocenters: (1R,5R,6R,7R,9S,11S,12S,13S,14S).
TETRODOTOXIN = (
    "C([C@@]1([C@H]2[C@@H]3[C@H](N=C(N[C@@]34[C@@H]([C@@H]1O[C@]"
    "([C@H]4O)(O2)O)O)N)O)O)O"
)

_BOND_ORDER = {
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.AROMATIC: 1.5,
}
_WEDGE = {Chem.BondDir.BEGINWEDGE: 1, Chem.BondDir.BEGINDASH: 2}


def _draw(smiles):
    """Lay out *smiles* in 2D with wedges, as a user would draw it."""
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    data = MolecularData()
    ids = []
    for atom in mol.GetAtoms():
        p = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        ids.append(
            data.add_atom(
                atom.GetSymbol(),
                (p.x / ANGSTROM_PER_PIXEL, -p.y / ANGSTROM_PER_PIXEL),
                atom.GetFormalCharge(),
            )
        )
    for bond in mol.GetBonds():
        data.add_bond(
            ids[bond.GetBeginAtomIdx()],
            ids[bond.GetEndAtomIdx()],
            _BOND_ORDER[bond.GetBondType()],
            _WEDGE.get(bond.GetBondDir(), 0),
        )
    return data


def _convert(data, mode="rdkit"):
    """Run the real 2D->3D pipeline and restore editor atom IDs, as the app does."""
    from moleditpy.ui.calculation_worker import CalculationWorker
    from moleditpy.ui.compute_logic import ComputeManager

    compute = ComputeManager.__new__(ComputeManager)
    compute.host = SimpleNamespace(state_manager=SimpleNamespace(data=data))
    block = ComputeManager._setup_mol_block_for_worker(
        compute, data.to_rdkit_mol(use_2d_stereo=False)
    )
    worker = CalculationWorker()
    results = []
    worker.finished.connect(results.append)
    worker.error.connect(lambda e: results.append(("error", e)))
    worker.run_calculation(
        block,
        {"conversion_mode": mode, "worker_id": 1, "optimization_method": "MMFF_RDKIT"},
    )
    _, mol = results[-1]
    assert isinstance(mol, Chem.Mol), results[-1]
    for idx, atom_id in enumerate(data.atoms):
        mol.GetAtomWithIdx(idx).SetIntProp("_original_atom_id", atom_id)
    return mol


def _mirror(mol):
    """Return *mol* reflected through the YZ plane: every center inverts."""
    mirrored = Chem.Mol(mol)
    conf = mirrored.GetConformer()
    for i in range(mirrored.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (-p.x, p.y, p.z))
    return mirrored


def test_tetrodotoxin_drawing_specifies_all_nine_centers():
    """Every stereocenter in PubChem's tetrodotoxin is read from the wedges."""
    assert len(drawn_chirality(_draw(TETRODOTOXIN))) == 9


def test_tetrodotoxin_conversion_keeps_every_center(app):
    """Converting tetrodotoxin to 3D keeps all nine drawn configurations."""
    data = _draw(TETRODOTOXIN)
    mol = _convert(data)

    assert find_chirality_mismatches(data, mol) == []
    assert actual_chirality(mol) == drawn_chirality(data)


def test_tetrodotoxin_mirrored_result_flags_every_center(app):
    """A 3D result with every center inverted is reported center by center."""
    data = _draw(TETRODOTOXIN)
    mismatches = find_chirality_mismatches(data, _mirror(_convert(data)))

    assert len(mismatches) == 9
    drawn = drawn_chirality(data)
    for m in mismatches:
        assert m.drawn == drawn[m.atom_id]
        assert m.actual == {"R": "S", "S": "R"}[m.drawn]
        assert m.symbol == "C"
        assert m.rdkit_index is not None


def test_single_inverted_center_is_the_only_one_reported():
    """Only the center whose configuration changed is reported."""
    data = _draw("C[C@@H](O)[C@H](N)C(=O)O")
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    # Read the wedges before adding H: AddHs puts new atoms at the origin,
    # which would scramble the geometry the wedges are read against.
    Chem.AssignChiralTypesFromBondDirs(mol)
    mol = Chem.AddHs(mol)
    first, _ = sorted(drawn_chirality(data))
    target = next(
        a for a in mol.GetAtoms() if a.GetIntProp("_original_atom_id") == first
    )
    target.InvertChirality()
    AllChem.EmbedMolecule(mol, randomSeed=11)

    mismatches = find_chirality_mismatches(data, mol)

    assert [m.atom_id for m in mismatches] == [first]


def test_unspecified_centers_are_not_compared():
    """A drawing without wedges specifies nothing, so nothing can mismatch."""
    data = _draw("CC(O)C(N)C(=O)O")
    mol = Chem.AddHs(data.to_rdkit_mol())
    AllChem.EmbedMolecule(mol, randomSeed=3)

    assert drawn_chirality(data) == {}
    assert find_chirality_mismatches(data, mol) == []


def test_lost_stereocenter_reported_without_3d_label():
    """A drawn center with no configuration in 3D is reported as actual None."""
    data = _draw("C[C@@H](O)CC")
    flat = Chem.AddHs(data.to_rdkit_mol())
    conf = Chem.Conformer(flat.GetNumAtoms())
    conf.Set3D(True)
    flat.AddConformer(conf, assignId=True)  # every atom at the origin

    mismatches = find_chirality_mismatches(data, flat)

    assert len(mismatches) == 1
    assert mismatches[0].actual is None


def test_structure_without_conformer_is_not_checked():
    """No coordinates means no 3D configuration to compare."""
    data = _draw("C[C@@H](O)CC")
    assert actual_chirality(data.to_rdkit_mol()) == {}


def test_rdkit_failure_is_logged_not_raised(caplog):
    """An RDKit error during the check is logged and treated as no mismatch."""
    from unittest.mock import patch

    data = _draw("C[C@@H](O)CC")
    with patch(
        "moleditpy.core.stereo_check.rdCIPLabeler.AssignCIPLabels",
        side_effect=RuntimeError("cip failed"),
    ):
        assert find_chirality_mismatches(data, Chem.Mol()) == []
    assert "Chirality check failed" in caplog.text


@pytest.mark.parametrize(
    "filename", sorted(os.listdir(os.path.dirname(core_pkg.__file__)))
)
def test_core_has_no_ui_dependency(filename):
    """core/ stays free of Qt, PyVista/VTK and the ui package."""
    if not filename.endswith(".py"):
        pytest.skip("not a module")
    path = os.path.join(os.path.dirname(core_pkg.__file__), filename)
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    forbidden = ("PyQt6", "PyQt5", "PySide6", "pyvista", "pyvistaqt", "vtk")
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [("." * node.level) + (node.module or "")]
        else:
            continue
        for name in names:
            assert name.lstrip(".").split(".")[0] not in forbidden, (
                f"{filename}: {name}"
            )
            assert not name.startswith("..ui"), f"{filename}: {name}"


# --- E/Z labels are CIP E/Z ------------------------------------------------


def _labelled_alkene(smiles, label):
    """Draw *smiles* without stereo and put E/Z label *label* (3=Z, 4=E) on C=C."""
    flat = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False)
    data = _draw(flat)
    double = next(k for k, b in data.bonds.items() if b["order"] == 2)
    data.bonds[double]["stereo"] = label
    return data, double


def _cip_ez(mol, atom_ids):
    """CIP E/Z the 3D structure has at the bond between two editor atom IDs."""
    from rdkit.Chem import rdCIPLabeler

    probe = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(probe)
    rdCIPLabeler.AssignCIPLabels(probe)
    for bond in probe.GetBonds():
        ids = {
            bond.GetBeginAtom().GetIntProp("_original_atom_id"),
            bond.GetEndAtom().GetIntProp("_original_atom_id"),
        }
        if ids == set(atom_ids) and bond.HasProp("_CIPCode"):
            return bond.GetProp("_CIPCode")
    return None


@pytest.mark.parametrize(
    "smiles",
    [
        "CC=CC",  # first heavy neighbor is the CIP-higher one
        "CC(Cl)=CC",  # methyl picked over Cl: used to come out inverted
        "ClC(Br)=C(C)CC",  # both ends need the CIP order, not the first neighbor
        "OC(C)=C(N)C",
    ],
)
@pytest.mark.parametrize("label,expected", [(3, "Z"), (4, "E")])
def test_ez_label_reaches_3d_as_cip_descriptor(app, smiles, label, expected):
    """A Z/E label means CIP Z/E in 3D, whichever neighbor is listed first."""
    data, double = _labelled_alkene(smiles, label)
    assert _cip_ez(_convert(data), double) == expected


def test_align_ez_to_cip_flips_only_wrong_bonds():
    """align_ez_to_cip flips a bond that reads the other way and leaves the rest."""
    from moleditpy.core.molecular_data import align_ez_to_cip
    from rdkit.Chem import rdCIPLabeler

    mol = Chem.MolFromSmiles("C/C(Cl)=C/C")  # trans methyls: CIP Z (Cl > CH3)
    bond = next(b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    before = bond.GetStereo()

    align_ez_to_cip(mol, {bond.GetIdx(): "Z"})
    assert bond.GetStereo() == before  # already Z: untouched

    align_ez_to_cip(mol, {bond.GetIdx(): "E"})
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    assert probe.GetBondWithIdx(bond.GetIdx()).GetProp("_CIPCode") == "E"


def test_align_ez_to_cip_ignores_bonds_without_cip_ez():
    """Two identical substituents on one end: no CIP E/Z, nothing to align."""
    from moleditpy.core.molecular_data import align_ez_to_cip

    mol = Chem.MolFromSmiles("CC(C)=CC")
    bond = next(b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    before = bond.GetStereo()
    align_ez_to_cip(mol, {bond.GetIdx(): "Z"})
    assert bond.GetStereo() == before


@pytest.mark.parametrize("stereo,expected", [("Z", "Z"), ("E", "E")])
def test_worker_explicit_stereo_uses_cip(stereo, expected):
    """The worker's M CFG fallback applies E/Z by CIP rank too."""
    from rdkit.Chem import AllChem, rdCIPLabeler

    from moleditpy.ui.calculation_worker import _apply_explicit_stereo

    mol = Chem.AddHs(Chem.MolFromSmiles("CC(Cl)=CC"))
    bond = next(b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    wanted = Chem.BondStereo.STEREOZ if stereo == "Z" else Chem.BondStereo.STEREOE
    _apply_explicit_stereo(mol, {bond.GetIdx(): wanted})
    AllChem.EmbedMolecule(mol, randomSeed=3)
    Chem.AssignStereochemistryFrom3D(mol)
    rdCIPLabeler.AssignCIPLabels(mol)
    assert mol.GetBondWithIdx(bond.GetIdx()).GetProp("_CIPCode") == expected


# --- E/Z check -------------------------------------------------------------


def test_ez_check_passes_a_correct_conversion(app):
    """A labelled double bond converted correctly reports nothing."""
    from moleditpy.core.stereo_check import find_ez_mismatches

    data, _ = _labelled_alkene("CC(Cl)=CC", 3)
    assert find_ez_mismatches(data, _convert(data)) == []


def test_ez_check_reports_a_wrong_double_bond(app):
    """A 3D bond with the other configuration is reported with its atoms."""
    from moleditpy.core.stereo_check import find_ez_mismatches

    data, double = _labelled_alkene("CC(Cl)=CC", 3)
    mol = _convert(data)
    data.bonds[double]["stereo"] = 4  # the label now asks for E; 3D is Z

    (m,) = find_ez_mismatches(data, mol)
    assert m.atom_ids == double
    assert (m.drawn, m.actual) == ("E", "Z")
    assert m.symbols == ("C", "C")
    assert m.rdkit_bond_index is not None
    assert m.rdkit_atom_indices is not None


def test_ez_labels_without_cip_ez_are_not_compared(app):
    """Two identical groups on one end: the label is unverifiable, skip it."""
    from moleditpy.core.stereo_check import drawn_ez, find_ez_mismatches

    data, _ = _labelled_alkene("CC(C)=CC", 3)
    assert drawn_ez(data) == {}
    assert find_ez_mismatches(data, _convert(data)) == []


def test_ez_check_needs_a_labelled_bond():
    """Without an E/Z label nothing is compared, even with a double bond."""
    from moleditpy.core.stereo_check import drawn_ez, find_ez_mismatches

    data = _draw("CC=CC")
    mol = Chem.AddHs(data.to_rdkit_mol())
    AllChem.EmbedMolecule(mol, randomSeed=2)
    assert drawn_ez(data) == {}
    assert find_ez_mismatches(data, mol) == []


def test_ez_check_without_conformer_is_skipped():
    """No 3D coordinates, nothing to compare."""
    from moleditpy.core.stereo_check import find_ez_mismatches

    data, _ = _labelled_alkene("CC(Cl)=CC", 3)
    assert find_ez_mismatches(data, data.to_rdkit_mol()) == []
