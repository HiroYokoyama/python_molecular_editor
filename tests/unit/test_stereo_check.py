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


def test_actual_ez_reads_cip_from_3d():
    """actual_ez labels each stereo double bond from the 3D coordinates."""
    from rdkit.Chem import AllChem
    from moleditpy.core.stereo_check import actual_ez

    for smiles, want in (("C/C=C/C", "E"), ("C/C=C\C", "Z")):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=7)
        double = next(
            b.GetIdx() for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 2
        )
        assert actual_ez(mol) == {double: want}

    flat = Chem.MolFromSmiles("C/C=C/C")  # no conformer
    assert actual_ez(flat) == {}


def test_drawn_stereo_is_empty_when_the_drawing_cannot_be_built(monkeypatch):
    """An unbuildable drawing specifies nothing to compare."""
    from moleditpy.core import stereo_check as sc

    data = _draw("C[C@H](F)Cl")
    monkeypatch.setattr(data, "to_rdkit_mol", lambda **_kw: None)
    assert sc.drawn_chirality(data) == {}

    alkene, _ = _labelled_alkene("CC=CC", 4)
    monkeypatch.setattr(alkene, "to_rdkit_mol", lambda **_kw: None)
    assert sc.drawn_ez(alkene) == {}


def test_3d_checks_skip_a_molecule_without_coordinates():
    """Without a conformer there is no 3D configuration to read."""
    from moleditpy.core import stereo_check as sc

    assert sc.actual_chirality(Chem.MolFromSmiles("C[C@H](F)Cl")) == {}
    data, _ = _labelled_alkene("CC=CC", 4)
    assert sc.find_ez_mismatches(data, Chem.MolFromSmiles("CC=CC")) == []


def test_ez_check_survives_a_cip_failure():
    """A CIP labeller error is logged and the check reports nothing."""
    from unittest.mock import patch

    from moleditpy.core import stereo_check as sc

    data, _ = _labelled_alkene("CC=CC", 4)
    mol = _convert(data)
    with patch.object(
        sc.rdCIPLabeler, "AssignCIPLabels", side_effect=RuntimeError("boom")
    ):
        assert sc.find_ez_mismatches(data, mol) == []


def test_ez_check_skips_labels_on_atoms_no_longer_drawn():
    """A label whose atoms were deleted since is not compared."""
    from unittest.mock import patch

    from moleditpy.core import stereo_check as sc

    data, _ = _labelled_alkene("CC=CC", 4)
    mol = _convert(data)
    with patch.object(sc, "drawn_ez", return_value={(998, 999): "E"}):
        assert sc.find_ez_mismatches(data, mol) == []


def _three_arm_center(cl_angle, c1_angle, c3_angle, stereo):
    """CH(Cl)(CH3)(C2H5) with Cl wedged (1) or hashed (2), arms at the given angles."""
    import math

    data = MolecularData()
    center = data.add_atom("C", (0.0, 0.0))

    def arm(symbol, angle):
        rad = math.radians(angle)
        return data.add_atom(symbol, (75 * math.cos(rad), -75 * math.sin(rad)))

    cl, c1, c3 = arm("Cl", cl_angle), arm("C", c1_angle), arm("C", c3_angle)
    x, y = data.atoms[c3]["pos"]
    c4 = data.add_atom("C", (x + 75, y))
    data.add_bond(center, cl, 1, stereo)
    data.add_bond(center, c1)
    data.add_bond(center, c3)
    data.add_bond(c3, c4)
    return data


@pytest.mark.parametrize("stereo", [1, 2])
@pytest.mark.parametrize("bend", [0.0, 2.0, -2.0])
def test_t_shaped_center_reads_like_the_bent_drawing(stereo, bend):
    """Plain bonds in a straight line read the same as the usual 120-degree layout.

    RDKit alone flags the straight (snapped) drawing as ambiguous or inverts
    it, while bending the plain bonds by two degrees reads correctly.
    """
    usual = drawn_chirality(_three_arm_center(90, 210, 330, stereo))
    t_shape = drawn_chirality(_three_arm_center(90, 180 + bend, 360 - bend, stereo))

    assert usual == {0: "R" if stereo == 1 else "S"}
    assert t_shape == usual


@pytest.mark.parametrize("stereo", [1, 2])
def test_t_shaped_center_converts_to_the_drawn_isomer(app, stereo):
    """The 3D structure from a T-shaped wedge drawing has the drawn configuration."""
    data = _three_arm_center(90, 180, 0, stereo)
    mol = _convert(data)

    assert find_chirality_mismatches(data, mol) == []
    assert actual_chirality(mol) == {0: "R" if stereo == 1 else "S"}


def test_wedge_inside_a_narrow_fan_is_left_to_rdkit():
    """A wedge between two plain bonds less than 180 degrees apart keeps RDKit's reading."""
    from moleditpy.core.molecular_data import wedge_reading

    data = _three_arm_center(90, 150, 30, 1)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    assert wedge_reading(mol, 0) is None
    probe = Chem.Mol(mol)
    Chem.AssignChiralTypesFromBondDirs(probe)
    assert drawn_chirality(data) == {0: "S"}


@pytest.mark.parametrize("smiles", [TETRODOTOXIN, "C[C@](F)(Cl)Br", "C[S@](=O)CC"])
def test_rdkit_wedged_drawings_are_left_unchanged(smiles):
    """Wedges RDKit places itself read the same before and after the fix."""
    from moleditpy.core.molecular_data import wedge_reading

    data = _draw(smiles)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    before = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondDir()) for b in mol.GetBonds()
    }
    probe = Chem.Mol(mol)
    Chem.AssignChiralTypesFromBondDirs(probe)
    for begin, _, direction in before:
        if direction in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            reading = wedge_reading(mol, begin)
            assert reading in (None, probe.GetAtomWithIdx(begin).GetChiralTag())
    assert {(i, j) for i, j, _ in before} == {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
    }


def _cross(down_stereo, symbols=("F", "Br", "Cl")):
    """C(CH3)(F)(Br)(Cl): methyl and F left/right, Br wedged up, Cl down."""
    data = MolecularData()
    center = data.add_atom("C", (0.0, 0.0))
    me = data.add_atom("C", (-75.0, 0.0))
    right = data.add_atom(symbols[0], (75.0, 0.0))
    up = data.add_atom(symbols[1], (0.0, -75.0))
    down = data.add_atom(symbols[2], (0.0, 75.0))
    data.add_bond(center, me)
    data.add_bond(center, right)
    data.add_bond(center, up, 1, 1)
    data.add_bond(center, down, 1, down_stereo)
    return data


def test_opposite_wedge_and_hash_is_reported_as_ambiguous():
    """A wedge and a hash on opposite bonds specify nothing and are flagged."""
    from moleditpy.core.stereo_check import ambiguous_wedge_centers

    data = _cross(2)
    assert drawn_chirality(data) == {}
    assert ambiguous_wedge_centers(data) == [0]


def test_single_wedge_cross_is_not_ambiguous():
    """One wedge on a four-arm cross is a readable drawing."""
    from moleditpy.core.stereo_check import ambiguous_wedge_centers

    data = _cross(0)
    assert drawn_chirality(data) == {0: "S"}
    assert ambiguous_wedge_centers(data) == []


def test_ambiguous_wedges_on_a_non_stereocenter_are_ignored():
    """Two identical substituents: no configuration to lose, nothing to report."""
    from moleditpy.core.stereo_check import ambiguous_wedge_centers

    assert ambiguous_wedge_centers(_cross(2, symbols=("C", "C", "Cl"))) == []


def test_find_ambiguous_centers_reports_the_3d_configuration(app):
    """The row carries "ambiguous" as drawn and whatever the 3D structure got."""
    from moleditpy.core.stereo_check import find_ambiguous_centers

    data = _cross(2)
    mol = _convert(data)
    rows = find_ambiguous_centers(data, mol)

    assert len(rows) == 1
    assert rows[0].atom_id == 0 and rows[0].drawn == "ambiguous"
    assert rows[0].actual in ("R", "S")
    assert find_ambiguous_centers(_cross(0), mol) == []


def test_wedge_reading_declines_unreadable_centers():
    """No conformer, wrong degree, overlapping atoms or no stereo bond: None."""
    from moleditpy.core.molecular_data import wedge_reading

    mol = _three_arm_center(90, 210, 330, 1).to_rdkit_mol(use_2d_stereo=False)
    assert wedge_reading(_no_conf(mol), 0) is None
    assert wedge_reading(mol, 1) is None  # Cl: one neighbor
    assert wedge_reading(mol, 3) is None  # ethyl CH2: no wedge starts here

    overlapping = Chem.Mol(mol)
    overlapping.GetConformer().SetAtomPosition(2, (0.0, 0.0, 0.0))
    assert wedge_reading(overlapping, 0) is None


def _no_conf(mol):
    copy = Chem.Mol(mol)
    copy.RemoveAllConformers()
    return copy


def test_resolve_wedge_chirality_skips_a_molecule_without_coordinates():
    """Without a conformer there is nothing to read; the wedges stay."""
    from moleditpy.core.molecular_data import resolve_wedge_chirality

    mol = _no_conf(_three_arm_center(90, 180, 0, 1).to_rdkit_mol(use_2d_stereo=False))
    before = [b.GetBondDir() for b in mol.GetBonds()]
    resolve_wedge_chirality(mol)
    assert [b.GetBondDir() for b in mol.GetBonds()] == before


def test_resolve_keeps_the_drawing_when_rewedging_does_not_hold(monkeypatch):
    """If RDKit cannot be made to read the intended center, the user's wedges stay."""
    import moleditpy.core.molecular_data as md

    data = _three_arm_center(90, 180, 0, 1)
    monkeypatch.setattr(md, "resolve_wedge_chirality", lambda mol: None)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    monkeypatch.undo()
    before = [(b.GetBeginAtomIdx(), b.GetBondDir()) for b in mol.GetBonds()]

    monkeypatch.setattr(md.Chem, "WedgeMolBonds", lambda mol, conf: None)
    md.resolve_wedge_chirality(mol)
    assert [(b.GetBeginAtomIdx(), b.GetBondDir()) for b in mol.GetBonds()] == before


def test_resolve_wedge_chirality_logs_rdkit_failure(monkeypatch, caplog):
    """An RDKit error while resolving is logged, not raised."""
    import moleditpy.core.molecular_data as md

    def boom(mol):
        raise RuntimeError("boom")

    mol = _three_arm_center(90, 180, 0, 1).to_rdkit_mol(use_2d_stereo=False)
    monkeypatch.setattr(md.Chem, "AssignChiralTypesFromBondDirs", boom)
    with caplog.at_level("WARNING"):
        md.resolve_wedge_chirality(mol)
    assert "Resolving wedge stereochemistry failed" in caplog.text


def test_ambiguous_check_skips_an_unbuildable_drawing(monkeypatch):
    """If the drawing cannot be built, nothing is reported."""
    from moleditpy.core.stereo_check import ambiguous_wedge_centers

    data = _cross(2)
    monkeypatch.setattr(data, "to_rdkit_mol", lambda use_2d_stereo=True: None)
    assert ambiguous_wedge_centers(data) == []


def test_ambiguous_check_survives_an_rdkit_failure(monkeypatch, caplog):
    """An RDKit error in the ambiguity check is logged and reports nothing."""
    import moleditpy.core.stereo_check as sc

    def boom(data):
        raise RuntimeError("boom")

    monkeypatch.setattr(sc, "ambiguous_wedge_centers", boom)
    with caplog.at_level("WARNING"):
        assert sc.find_ambiguous_centers(_cross(2), Chem.Mol()) == []
    assert "Ambiguous wedge check failed" in caplog.text


def test_wedge_reading_needs_a_net_wedge_on_three_neighbors():
    """No stereo bond, or a wedge and a hash that cancel, leave the center unread."""
    from moleditpy.core.molecular_data import wedge_reading

    plain = _three_arm_center(90, 210, 330, 0).to_rdkit_mol(use_2d_stereo=False)
    assert wedge_reading(plain, 0) is None

    data = _three_arm_center(90, 210, 330, 1)
    data.add_bond(0, 2, 1, 2)
    assert wedge_reading(data.to_rdkit_mol(use_2d_stereo=False), 0) is None


def test_resolve_rejects_rewedging_that_touches_other_bonds(monkeypatch):
    """A re-wedge that alters a bond away from the fixed center is not applied."""
    import moleditpy.core.molecular_data as md

    data = _three_arm_center(90, 180, 0, 1)
    monkeypatch.setattr(md, "resolve_wedge_chirality", lambda mol: None)
    mol = data.to_rdkit_mol(use_2d_stereo=False)
    monkeypatch.undo()
    before = [(b.GetBeginAtomIdx(), b.GetBondDir()) for b in mol.GetBonds()]

    real_wedge = md.Chem.WedgeMolBonds

    def wedge_and_stray(target, conf):
        real_wedge(target, conf)
        # C3-C4 does not touch the T-shaped center (atom 0).
        target.GetBondBetweenAtoms(3, 4).SetBondDir(md.Chem.BondDir.BEGINWEDGE)

    monkeypatch.setattr(md.Chem, "WedgeMolBonds", wedge_and_stray)
    md.resolve_wedge_chirality(mol)
    assert [(b.GetBeginAtomIdx(), b.GetBondDir()) for b in mol.GetBonds()] == before
