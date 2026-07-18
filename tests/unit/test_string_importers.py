"""Tests for SMILES/InChI import — verifying atom counts, bonds, and properties against RDKit reference."""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from moleditpy.ui.string_importers import StringImporterManager
from PyQt6.QtCore import QPointF, QTimer
from unittest.mock import patch


class DummyImporter(StringImporterManager):
    """Thin wrapper to call importer methods with a mocked host."""

    def __init__(self, host):
        self._host = host
        self.host = host  # required by StringImporterManager (manager architecture)
        self.data = host.state_manager.data
        self.scene = host.init_manager.scene

    def __getattr__(self, name):
        return getattr(self._host, name)

    def check_unsaved_changes(self):
        return True  # Always proceed

    def clear_2d_editor(self, push_to_undo=True):
        pass

    def statusBar(self):
        return self._host.statusBar()


# =============================================================================
# SMILES import tests — compare scene atoms against RDKit reference
# =============================================================================


def _load_and_get_data(mock_parser_host, importer_method, input_str):
    """Helper: run an importer, return the final MolecularData state."""
    importer = DummyImporter(mock_parser_host)
    importer.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    importer.view_2d.mapToScene.return_value = QPointF(0, 0)

    with patch.object(QTimer, "singleShot"):
        getattr(importer, importer_method)(input_str)

    return importer.data


def test_smiles_ethanol_atom_count(mock_parser_host):
    """Ethanol SMILES should produce 3 heavy atoms matching RDKit reference."""
    data = _load_and_get_data(mock_parser_host, "load_from_smiles", "CCO")

    ref_mol = Chem.MolFromSmiles("CCO")
    assert len(data.atoms) == ref_mol.GetNumAtoms()


def test_smiles_ethanol_bond_count(mock_parser_host):
    """Ethanol SMILES should produce correct bond count."""
    data = _load_and_get_data(mock_parser_host, "load_from_smiles", "CCO")

    ref_mol = Chem.MolFromSmiles("CCO")
    assert len(data.bonds) == ref_mol.GetNumBonds()


def test_smiles_benzene(mock_parser_host):
    """Benzene SMILES should produce 6C atoms and 6 bonds."""
    data = _load_and_get_data(mock_parser_host, "load_from_smiles", "c1ccccc1")

    assert len(data.atoms) == 6
    assert len(data.bonds) == 6
    assert all(a["symbol"] == "C" for a in data.atoms.values())


def test_smiles_preserves_formal_charge(mock_parser_host):
    """Ammonium NH4+ should carry charge +1 on nitrogen."""
    data = _load_and_get_data(mock_parser_host, "load_from_smiles", "[NH4+]")

    n_atoms = [a for a in data.atoms.values() if a["symbol"] == "N"]
    assert len(n_atoms) == 1
    assert n_atoms[0]["charge"] == 1


def test_smiles_aspirin_formula(mock_parser_host):
    """Aspirin SMILES should produce correct molecular formula vs RDKit."""
    data = _load_and_get_data(
        mock_parser_host, "load_from_smiles", "CC(=O)Oc1ccccc1C(=O)O"
    )

    # Build RDKit mol from the imported data for formula comparison
    mol = data.to_rdkit_mol()
    ref = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)


def test_smiles_stereocenter_wedge(mock_parser_host):
    """Chiral SMILES should produce at least one wedge or dash bond."""
    data = _load_and_get_data(mock_parser_host, "load_from_smiles", "[C@@H](F)(Cl)Br")

    stereo_bonds = [b for b in data.bonds.values() if b.get("stereo", 0) in [1, 2]]
    assert len(stereo_bonds) >= 1


def test_smiles_invalid_shows_error(mock_parser_host):
    """Invalid SMILES should trigger status bar error, not crash."""
    importer = DummyImporter(mock_parser_host)
    importer.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    importer.view_2d.mapToScene.return_value = QPointF(0, 0)

    with patch.object(QTimer, "singleShot"):
        importer.load_from_smiles("INVALID_SMILES_XYZ")

    # Should call statusBar with specific SMILES error prefix
    mock_parser_host.statusBar().showMessage.assert_called()
    last_msg = mock_parser_host.statusBar().showMessage.call_args[0][0]
    assert last_msg.startswith("Invalid SMILES:"), (
        f"Expected 'Invalid SMILES:' prefix, got: {last_msg}"
    )


def test_smiles_empty_shows_error(mock_parser_host):
    """Empty SMILES should trigger error message."""
    importer = DummyImporter(mock_parser_host)
    with patch.object(QTimer, "singleShot"):
        importer.load_from_smiles("")

    mock_parser_host.statusBar().showMessage.assert_called()


# =============================================================================
# InChI import tests
# =============================================================================


def test_inchi_ethanol(mock_parser_host):
    """InChI import should match RDKit reference atom/bond count."""
    inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    data = _load_and_get_data(mock_parser_host, "load_from_inchi", inchi)

    ref = Chem.MolFromInchi(inchi)
    assert len(data.atoms) == ref.GetNumAtoms()
    assert len(data.bonds) == ref.GetNumBonds()


def test_inchi_caffeine_formula(mock_parser_host):
    """Caffeine InChI should produce correct molecular formula."""
    inchi = "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"
    data = _load_and_get_data(mock_parser_host, "load_from_inchi", inchi)

    mol = data.to_rdkit_mol()
    ref = Chem.MolFromInchi(inchi)
    assert rdMolDescriptors.CalcMolFormula(mol) == rdMolDescriptors.CalcMolFormula(ref)


def test_inchi_invalid_shows_error(mock_parser_host):
    """Invalid InChI should show error, not crash."""
    importer = DummyImporter(mock_parser_host)
    with patch.object(QTimer, "singleShot"):
        importer.load_from_inchi("InChI=INVALID")

    mock_parser_host.statusBar().showMessage.assert_called()
    last_msg = mock_parser_host.statusBar().showMessage.call_args[0][0]
    assert last_msg.startswith("Invalid InChI:"), (
        f"Expected 'Invalid InChI:' prefix, got: {last_msg}"
    )


# =============================================================================
# Input dialogs, placement, stereo mapping, and error paths
# =============================================================================

from unittest.mock import MagicMock


def test_smiles_dialog_ok_loads(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.load_from_smiles = MagicMock()
    with patch(
        "moleditpy.ui.string_importers.QInputDialog.getText",
        return_value=("CCO", True),
    ):
        importer.import_smiles_dialog()
    importer.load_from_smiles.assert_called_once_with("CCO")


def test_smiles_dialog_cancel_does_not_load(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.load_from_smiles = MagicMock()
    with patch(
        "moleditpy.ui.string_importers.QInputDialog.getText",
        return_value=("CCO", False),
    ):
        importer.import_smiles_dialog()
    importer.load_from_smiles.assert_not_called()


def test_smiles_dialog_empty_text_does_not_load(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.load_from_smiles = MagicMock()
    with patch(
        "moleditpy.ui.string_importers.QInputDialog.getText",
        return_value=("", True),
    ):
        importer.import_smiles_dialog()
    importer.load_from_smiles.assert_not_called()


def test_inchi_dialog_ok_loads(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.load_from_inchi = MagicMock()
    with patch(
        "moleditpy.ui.string_importers.QInputDialog.getText",
        return_value=("InChI=1S/CH4/h1H4", True),
    ):
        importer.import_inchi_dialog()
    importer.load_from_inchi.assert_called_once_with("InChI=1S/CH4/h1H4")


def test_inchi_dialog_cancel_does_not_load(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.load_from_inchi = MagicMock()
    with patch(
        "moleditpy.ui.string_importers.QInputDialog.getText",
        return_value=("x", False),
    ):
        importer.import_inchi_dialog()
    importer.load_from_inchi.assert_not_called()


def test_placement_center_offsets_right_of_existing_atoms(mock_parser_host):
    """Second import must land 80px right of the rightmost atom at mean y."""
    importer = DummyImporter(mock_parser_host)
    importer.data.atoms = {
        1: {"pos": QPointF(10.0, 20.0)},
        2: {"pos": QPointF(110.0, 40.0)},
    }
    center = importer._placement_center()
    assert center.x() == 190.0  # 110 + 80
    assert center.y() == 30.0  # mean of 20 and 40


def test_placement_center_supports_tuple_positions(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.data.atoms = {1: {"pos": (50.0, 10.0)}}
    center = importer._placement_center()
    assert center.x() == 130.0
    assert center.y() == 10.0


def test_placement_center_uses_viewport_when_empty(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer.data.atoms = {}
    importer.view_2d.mapToScene.return_value = QPointF(7.0, 8.0)
    center = importer._placement_center()
    assert center == QPointF(7.0, 8.0)


def _double_bond_stereo_values(mock_parser_host, smiles):
    importer = DummyImporter(mock_parser_host)
    importer.view_2d.viewport().rect().center.return_value = QPointF(0, 0)
    importer.view_2d.mapToScene.return_value = QPointF(0, 0)
    with patch.object(QTimer, "singleShot"):
        importer.load_from_smiles(smiles)
    return [
        c.kwargs.get("bond_stereo")
        for c in mock_parser_host.init_manager.scene.create_bond.call_args_list
    ]


def test_smiles_e_double_bond_maps_stereo_4(mock_parser_host):
    stereos = _double_bond_stereo_values(mock_parser_host, "C/C=C/C")
    assert 4 in stereos


def test_smiles_z_double_bond_maps_stereo_3(mock_parser_host):
    stereos = _double_bond_stereo_values(mock_parser_host, r"C/C=C\C")
    assert 3 in stereos


def _last_status(mock_parser_host):
    return str(mock_parser_host.statusBar().showMessage.call_args.args[0])


def test_smiles_truly_empty_after_strip_reports_empty(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    with patch("moleditpy.ui.string_importers.Chem.MolFromSmiles", return_value=None):
        importer.load_from_smiles("   ")
    assert "SMILES string was empty." in _last_status(mock_parser_host)


def test_smiles_parse_runtime_error_reported(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    with patch(
        "moleditpy.ui.string_importers.AllChem.Compute2DCoords",
        side_effect=RuntimeError("rdkit exploded"),
    ):
        importer.load_from_smiles("CCO")
    assert "Error parsing SMILES: rdkit exploded" in _last_status(mock_parser_host)


def test_smiles_placement_error_reported(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer._placement_center = MagicMock(side_effect=RuntimeError("view gone"))
    with patch.object(QTimer, "singleShot"):
        importer.load_from_smiles("CCO")
    assert "Error loading from SMILES: view gone" in _last_status(mock_parser_host)


def test_inchi_truly_empty_after_strip_reports_empty(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    with patch("moleditpy.ui.string_importers.Chem.MolFromInchi", return_value=None):
        importer.load_from_inchi("   ")
    assert "InChI string was empty." in _last_status(mock_parser_host)


def test_inchi_parse_runtime_error_reported(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    with patch(
        "moleditpy.ui.string_importers.AllChem.Compute2DCoords",
        side_effect=RuntimeError("rdkit exploded"),
    ):
        importer.load_from_inchi("InChI=1S/CH4/h1H4")
    assert "Error parsing InChI: rdkit exploded" in _last_status(mock_parser_host)


def test_inchi_placement_error_reported(mock_parser_host):
    importer = DummyImporter(mock_parser_host)
    importer._placement_center = MagicMock(side_effect=RuntimeError("view gone"))
    with patch.object(QTimer, "singleShot"):
        importer.load_from_inchi("InChI=1S/CH4/h1H4")
    assert "Error loading from InChI: view gone" in _last_status(mock_parser_host)
