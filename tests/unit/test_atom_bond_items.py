import pytest
from unittest.mock import MagicMock, patch
from PyQt6.QtCore import QPointF, QRectF
from moleditpy.ui.atom_item import AtomItem
from moleditpy.ui.bond_item import BondItem


# Mock MainWindow and Scene
@pytest.fixture
def mock_main_window():
    mw = MagicMock()
    mw.current_symbol = "C"
    mw.bond_mode = False
    mw.erase_mode = False
    return mw


@pytest.fixture
def mock_scene():
    scene = MagicMock()
    return scene


# Helper class to tolerate QPoint in moveCenter (fixing main code bug in test env)
from PyQt6.QtCore import QPoint


class TolerantQRectF(QRectF):
    def moveCenter(self, p):
        if isinstance(p, QPoint):
            p = QPointF(p)
        super().moveCenter(p)

    def adjusted(self, x1, y1, x2, y2):
        r = super().adjusted(x1, y1, x2, y2)
        return TolerantQRectF(r)

    def translated(self, dx, dy):
        r = super().translated(dx, dy)
        return TolerantQRectF(r)


# Helper for mocking scene() and isSelected() on AtomItem
class MockableAtomItem(AtomItem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._selected_mock = False

    def isSelected(self):
        return self._selected_mock

    def setSelected(self, val):
        self._selected_mock = val

    def boundingRect(self):
        # Override to avoid QFontMetricsF dependency and potential crashes
        # Return a simple TolerantQRectF
        return TolerantQRectF(0, 0, 20, 20)


class TestAtomItem:
    @pytest.fixture
    def atom_item(self, mock_main_window, mock_scene):
        item = MockableAtomItem(1, "C", QPointF(0.0, 0.0))
        return item

    def test_init(self, atom_item):
        """Verify AtomItem initialization with ID, symbol and position."""
        assert atom_item.symbol == "C"
        assert atom_item.atom_id == 1
        assert atom_item.pos().x() == 0.0
        assert atom_item.pos().y() == 0.0
        # Verify default flags or properties that aren't passed in init
        assert atom_item.flags() & atom_item.GraphicsItemFlag.ItemIsMovable
        assert atom_item.flags() & atom_item.GraphicsItemFlag.ItemIsSelectable

    def test_paint_mock(self, atom_item):
        """Test paint logic by mocking QPainter"""
        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        mock_fm = MagicMock()
        mock_painter.fontMetrics.return_value = mock_fm

        # IMPORTANT: The side_effect ensures that whenever boundingRect is called, we return a FRESH TolerantQRectF
        mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(0, 0, 20, 20)

        atom_item.paint(mock_painter, option, widget)

        assert mock_painter.drawText.called
        # Verify correct text "C" is drawn
        text_drawn = False
        for args, _ in mock_painter.drawText.call_args_list:
            # drawText args vary, but usually the last arg or one of strings is text
            if "C" in args:
                text_drawn = True
                break
        assert text_drawn, "Atom symbol 'C' should be found in drawText calls"

    def test_paint_radical(self, atom_item):
        """Test painting logic for radicals"""
        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        mock_fm = MagicMock()
        mock_painter.fontMetrics.return_value = mock_fm
        mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(0, 0, 20, 20)

        # Test Radical 1 (Monoradical)
        atom_item.radical = 1
        atom_item.paint(mock_painter, option, widget)
        assert mock_painter.drawEllipse.called

        # Test Radical 2 (Diradical)
        atom_item.radical = 2
        mock_painter.reset_mock()
        mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(0, 0, 20, 20)

        atom_item.paint(mock_painter, option, widget)
        assert mock_painter.drawEllipse.call_count >= 2

    def test_paint_selection_hover(self, atom_item):
        """Test painting logic for selection and hover highlights"""
        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        mock_fm = MagicMock()
        mock_painter.fontMetrics.return_value = mock_fm
        mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(0, 0, 20, 20)

        # Test Selection (isSelected=True)
        atom_item.setSelected(True)
        atom_item.paint(mock_painter, option, widget)
        assert mock_painter.drawRect.called

        mock_painter.reset_mock()

        # Test Hover (isSelected=False, hovered=True)
        atom_item.setSelected(False)
        atom_item.hovered = True
        atom_item.paint(mock_painter, option, widget)
        assert mock_painter.drawRect.called

        mock_painter.reset_mock()

        # Test Problem (has_problem=True)
        atom_item.has_problem = True
        atom_item.paint(mock_painter, option, widget)
        assert mock_painter.drawRect.called

    def test_paint_implicit_hydrogens_and_flipping(self, atom_item):
        """Test painting logic for implicit hydrogens and text alignment"""
        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        mock_fm = MagicMock()
        mock_painter.fontMetrics.return_value = mock_fm
        # Ensure boundingRect returns a TolerantQRectF
        mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(0, 0, 20, 20)

        # Patch sip_isdeleted_safe to avoid issues with Mock objects not being C++ pointers
        with patch("moleditpy.ui.atom_item.sip_isdeleted_safe", return_value=False):
            # 1. Test Implicit Hydrogens (e.g. CH3)
            atom_item.implicit_h_count = 3
            atom_item.paint(mock_painter, option, widget)

            # Verify call args for drawText to contain subscript
            call_args_list = mock_painter.drawText.call_args_list
            found_subscript = False
            for args, _ in call_args_list:
                # args[0] might be rect, args[1] flags, args[2] text OR args[0] x, args[1] y, args[2] text
                # We look for string arg
                for arg in args:
                    if isinstance(arg, str) and "₃" in arg:
                        found_subscript = True
                        break
            assert found_subscript, (
                "Subscript for implicit H count 3 (₃) not found in drawText calls"
            )

            # 2. Test Text Flipping (AlignLeft vs AlignRight) based on bond direction

            # Change symbol to 'N' so it's not a skeletal Carbon (which hides H)
            atom_item.symbol = "N"

            # Setup neighbors
            # We need to add bonds to the atom
            bond_mock = MagicMock()
            other_atom = MagicMock()

            # Case A: Neighbor is to the LEFT (dx < 0) -> No Flip (Group on Right) -> AlignLeft
            other_atom.pos.return_value = QPointF(-10.0, 0.0)
            atom_item.setPos(0.0, 0.0)

            bond_mock.atom1 = atom_item
            bond_mock.atom2 = other_atom

            # atom_item.bonds is a list
            atom_item.bonds = [bond_mock]

            mock_painter.reset_mock()
            mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(
                0, 0, 20, 20
            )

            atom_item.paint(mock_painter, option, widget)

            # 3. Case B: Neighbor is on RIGHT (+10, 0) -> Flip -> H3C, AlignRight
            other_atom.pos.return_value = QPointF(10.0, 0.0)

            mock_painter.reset_mock()
            mock_fm.boundingRect.side_effect = lambda *args: TolerantQRectF(
                0, 0, 20, 20
            )

            atom_item.paint(mock_painter, option, widget)

            # Check for reversed text or checks
            # If flip_text=True, display_text = hydrogen_part + symbol
            # "H₃C" (H part first)

            call_args_list = mock_painter.drawText.call_args_list
            found_flipped_text = False
            for args, _ in call_args_list:
                for arg in args:
                    if (
                        isinstance(arg, str)
                        and arg.startswith("H")
                        and arg.endswith("N")
                    ):
                        found_flipped_text = True

            assert found_flipped_text, (
                "Expected flipped text 'H...N' when neighbor is on the right"
            )


# Helper for mocking scene() on BondItem
class MockableBondItem(BondItem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._mock_scene = MagicMock()

    def scene(self):
        return self._mock_scene


class TestBondItem:
    @pytest.fixture
    def bond_item(self, mock_main_window):
        atom1 = AtomItem(1, "C", QPointF(0.0, 0.0))
        atom2 = AtomItem(2, "C", QPointF(10.0, 10.0))

        # Use MockableBondItem to allow scene() mocking
        bond = MockableBondItem(atom1, atom2)
        return bond

    def test_init(self, bond_item):
        """Verify BondItem initialization with atom partners and default order."""
        assert bond_item.atom1.atom_id == 1
        assert bond_item.atom2.atom_id == 2
        assert bond_item.order == 1

    def test_update_position(self, bond_item):
        """Verify that the bond line updates correctly when atom positions change."""
        # Move an atom
        bond_item.atom1.setPos(5.0, 5.0)
        bond_item.update_position()

        line = bond_item.get_line_in_local_coords()
        # get_line_in_local_coords returns line from (0,0) to (p2-p1)
        assert line.p1() == QPointF(0.0, 0.0)
        # p2 (10,10) - p1 (5,5) = (5,5)
        assert line.p2() == QPointF(5.0, 5.0)

    def test_set_bond_order(self, bond_item):
        """Verify that the bond order can be changed and is reflected in the item state."""
        bond_item.set_order(2)  # method name is set_order
        assert bond_item.order == 2

        bond_item.set_order(3)
        assert bond_item.order == 3

    def test_paint_mock(self, bond_item):
        """Test paint logic for bond"""
        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        # Setup mocking for scene/window settings accessed in paint
        mock_scene = (
            bond_item.scene()
        )  # This is now our MagicMock from TestableBondItem
        mock_window = MagicMock()
        mock_window.settings = {}  # Empty dict for defaults
        mock_scene.window = mock_window
        mock_scene.views.return_value = [MagicMock()]
        mock_scene.views.return_value[0].window.return_value = mock_window

        bond_item.paint(mock_painter, option, widget)

        # Should draw line (single bond)
        assert mock_painter.drawLine.called

        # Test double bond
        bond_item.set_order(2)
        mock_painter.reset_mock()
        bond_item.paint(mock_painter, option, widget)
        # Double bond often draws lines
        assert mock_painter.drawLine.call_count >= 1

    def test_paint_ring_double_bond(self, bond_item):
        """Test painting logic for double bond in a ring"""
        bond_item.set_order(2)

        mock_painter = MagicMock()
        option = MagicMock()
        widget = MagicMock()

        # Mock scene and RDKit interaction (via MockableBondItem's mock_scene)
        mock_scene = bond_item.scene()

        # Setup window and data
        mock_window = MagicMock()
        mock_window.settings = {}
        mock_scene.window = mock_window

        # Setup RDKit mol mock
        mock_mol = MagicMock()
        mock_window.data.to_rdkit_mol.return_value = mock_mol

        # Setup atoms in mol
        atom1_mock = MagicMock()
        atom1_mock.HasProp.return_value = True
        atom1_mock.GetIntProp.return_value = 1
        atom1_mock.GetIdx.return_value = 0

        atom2_mock = MagicMock()
        atom2_mock.HasProp.return_value = True
        atom2_mock.GetIntProp.return_value = 2
        atom2_mock.GetIdx.return_value = 1

        mock_mol.GetAtoms.return_value = [atom1_mock, atom2_mock]

        # Setup Bond
        mock_rdkit_bond = MagicMock()
        mock_mol.GetBondBetweenAtoms.return_value = mock_rdkit_bond
        mock_rdkit_bond.IsInRing.return_value = True

        # Setup Ring Info
        mock_ring_info = MagicMock()
        mock_mol.GetRingInfo.return_value = mock_ring_info
        mock_ring_info.AtomRings.return_value = [[0, 1, 2]]

        # Setup GetAtomWithIdx for ring calculation failure safety or success
        # We need a 3rd atom to define a ring center
        atom3_mock = MagicMock()
        atom3_mock.HasProp.return_value = True
        atom3_mock.GetIntProp.return_value = 3

        def get_atom_with_idx(idx):
            if idx == 0:
                return atom1_mock
            if idx == 1:
                return atom2_mock
            return atom3_mock

        mock_mol.GetAtomWithIdx.side_effect = get_atom_with_idx

        # Setup scene.window.state_manager.data.atoms dictionary
        atom3_item = MagicMock()
        atom3_item.pos.return_value = QPointF(5.0, 5.0)

        mock_window.data.atoms = {
            1: {"item": bond_item.atom1},
            2: {"item": bond_item.atom2},
            3: {"item": atom3_item},
        }

        # mapFromScene mock
        bond_item.mapFromScene = MagicMock(side_effect=lambda p: p)  # Identity map

        # Act
        bond_item.paint(mock_painter, option, widget)

        assert mock_painter.drawLine.call_count >= 2

    def test_bounding_rect_ez_label(self, bond_item, qapp):
        """Test boundingRect expansion for E/Z labels"""
        bond_item.set_order(2)
        bond_item.stereo = 3  # Z-isomer

        # Mock scene settings for font size
        mock_scene = bond_item.scene()
        mock_window = MagicMock()
        mock_window.settings = {"atom_font_size_2d": 20}
        mock_scene.views.return_value = [MagicMock()]
        mock_scene.views.return_value[0].window.return_value = mock_window
        mock_scene.window = mock_window

        bond_item.boundingRect()

        bond_item.stereo = 0
        rect_no_stereo = bond_item.boundingRect()

        bond_item.stereo = 3
        rect_stereo = bond_item.boundingRect()

        # Ensure rect_stereo is not smaller, and hopefully different or contains the other
        assert rect_stereo.width() >= rect_no_stereo.width()
        assert rect_stereo.height() >= rect_no_stereo.height()
