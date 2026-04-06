import pytest
from unittest.mock import MagicMock, Mock, patch, mock_open
import numpy as np
from PyQt6.QtWidgets import QMainWindow
from moleditpy.ui.export_logic import ExportManager
from PyQt6.QtCore import QRectF


# Create dummy classes for reliable isinstance checks
class DummyAtomItem:
    def isVisible(self):
        return True

    def sceneBoundingRect(self):
        return QRectF(0, 0, 100, 100)

    def hide(self):
        pass

    def setVisible(self, v):
        pass


class DummyBondItem:
    def isVisible(self):
        return True

    def sceneBoundingRect(self):
        return QRectF(0, 0, 100, 100)

    def hide(self):
        pass

    def setVisible(self, v):
        pass


class TinyMesh:
    """Very simple class to avoid mock equality/isinstance issues."""
    def __init__(self, n_points=0):
        self.n_points = n_points
        self.point_data = {}
    def copy(self):
        # Return a NEW TinyMesh with same points
        return TinyMesh(self.n_points)
    def merge(self, other):
        return TinyMesh(self.n_points + getattr(other, "n_points", 0))
    def save(self, *args, **kwargs):
        pass
    def extract_points(self, *args, **kwargs):
        return TinyMesh(self.n_points // 2)


# Create a mock class that uses ExportManager
class MockMainWindow(ExportManager, QMainWindow):
    def __init__(self):
        # Initialize ExportManager with self as host
        ExportManager.__init__(self, self)

        # Initialize all managers to avoid AttributeErrors
        self.view_3d_manager = MagicMock()
        self.edit_3d_manager = MagicMock()
        self.state_manager = MagicMock()
        self.init_manager = MagicMock()
        self.io_manager = MagicMock()
        self.compute_manager = MagicMock()
        self.ui_manager = MagicMock()
        self.export_manager = self

        self.view_3d_manager.current_mol = MagicMock()
        self.host.init_manager.current_file_path = "/path/to/molecule.xyz"
        self.view_3d_manager.plotter = MagicMock()
        self.init_manager.scene = MagicMock()
        self._status_bar = MagicMock()
        self.init_manager.view_2d = MagicMock()
        self.state_manager.data = MagicMock()

    def statusBar(self):
        return self._status_bar


@pytest.fixture
def window():
    return MockMainWindow()


@pytest.fixture
def mock_file_dialog():
    with patch("moleditpy.ui.export_logic.QFileDialog") as mock:
        yield mock


@pytest.fixture
def mock_message_box():
    with patch("moleditpy.ui.export_logic.QMessageBox") as mock:
        yield mock


@pytest.fixture
def mock_pv_polydata():
    with patch("moleditpy.ui.export_logic.pv.PolyData", create=True) as mock:
        yield mock


def test_export_stl_success(window, mock_file_dialog, mock_pv_polydata):
    """Test export_stl success path."""
    mock_file_dialog.getSaveFileName.return_value = (
        "/path/to/export.stl",
        "STL Files (*.stl)",
    )

    mock_combined_mesh = MagicMock()
    mock_combined_mesh.n_points = 100

    with patch.object(
        window, "export_from_3d_view_no_color", return_value=mock_combined_mesh
    ):
        window.export_stl()

    mock_combined_mesh.save.assert_called_once_with("/path/to/export.stl", binary=True)
    window.statusBar().showMessage.assert_any_call(
        "STL exported to /path/to/export.stl"
    )


def test_export_stl_cancel(window, mock_file_dialog):
    """Test export_stl cancellation."""
    mock_file_dialog.getSaveFileName.return_value = ("", "")
    window.export_stl()
    call_args_list = window.statusBar().showMessage.call_args_list
    assert not any("STL exported" in str(args) for args, _ in call_args_list)


def test_export_stl_no_molecule(window):
    """Test export_stl with no molecule."""
    window.view_3d_manager.current_mol = None
    window.export_stl()
    window.statusBar().showMessage.assert_called_with(
        "Error: Please generate a 3D structure first."
    )


def test_export_obj_mtl_success(window, mock_file_dialog):
    """Test export_obj_mtl success path."""
    mock_file_dialog.getSaveFileName.return_value = (
        "/path/to/export.obj",
        "OBJ Files (*.obj)",
    )

    meshes = [
        {
            "mesh": MagicMock(),
            "color": [255, 0, 0],
            "name": "test_mesh",
            "type": "actor",
            "actor_name": "test",
        }
    ]

    with (
        patch.object(window, "export_from_3d_view_with_colors", return_value=meshes),
        patch.object(window, "create_multi_material_obj") as mock_create,
    ):
        window.export_obj_mtl()

        mock_create.assert_called_once_with(
            meshes, "/path/to/export.obj", "/path/to/export.mtl"
        )
        window.statusBar().showMessage.assert_any_call(
            "OBJ+MTL files with individual colors exported to /path/to/export.obj and /path/to/export.mtl"
        )


def test_export_color_stl_success(window, mock_file_dialog):
    """Test export_color_stl success path."""
    mock_file_dialog.getSaveFileName.return_value = (
        "/path/to/color_export.stl",
        "STL Files (*.stl)",
    )

    mock_combined_mesh = MagicMock()
    mock_combined_mesh.n_points = 100

    with patch.object(window, "export_from_3d_view", return_value=mock_combined_mesh):
        window.export_color_stl()

    mock_combined_mesh.save.assert_called_once_with(
        "/path/to/color_export.stl", binary=True
    )
    window.statusBar().showMessage.assert_any_call(
        "STL exported to /path/to/color_export.stl"
    )


def test_create_multi_material_obj_logic(window):
    """Test the file writing logic of create_multi_material_obj."""
    mock_mesh = MagicMock()
    mock_mesh.n_points = 3
    mock_mesh.points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    mock_mesh.n_cells = 1
    mock_cell = MagicMock()
    mock_cell.type = 5  # VTK_TRIANGLE
    mock_cell.point_ids = [0, 1, 2]
    mock_mesh.get_cell.return_value = mock_cell

    meshes = [
        {
            "mesh": mock_mesh,
            "color": [255, 0, 0],
            "name": "test_mesh",
            "type": "actor",
            "actor_name": "test",
        }
    ]

    obj_path = "test.obj"
    mtl_path = "test.mtl"

    with patch("builtins.open", mock_open()) as mock_file:
        window.create_multi_material_obj(meshes, obj_path, mtl_path)
        mock_file.assert_any_call(obj_path, "w")
        mock_file.assert_any_call(mtl_path, "w")

        handle = mock_file()
        write_calls = handle.write.call_args_list
        content = "".join([args[0] for args, _ in write_calls])

        assert "mtllib test.mtl" in content
        assert "usemtl material_0_test_mesh" in content
        assert "v 0.000000 0.000000 0.000000" in content
        assert "f 1 2 3" in content


def test_export_from_3d_view_logic(window):
    """Test export_from_3d_view to ensure it iterates actors and extracts meshes."""
    mock_mesh_data = TinyMesh(10)

    # actor setup using Mock to avoid auto-creating methods like extract_surface
    mock_actor = Mock(name="Actor")
    mock_actor.mapper = Mock()
    mock_actor.mapper.input = mock_mesh_data
    mock_actor.prop = Mock()
    mock_actor.prop.color = [1.0, 1.0, 1.0]

    window.view_3d_manager.plotter.renderer.actors = {"actor1": mock_actor}
    window.view_3d_manager.plotter.mesh = {"actor1": mock_mesh_data}

    with patch("moleditpy.ui.export_logic.pv") as mock_pv:
        mock_pv.PolyData = TinyMesh # Class-like
        mock_pv.wrap.side_effect = lambda x: x # Wrap returns input if already TinyMesh

        result = window.export_from_3d_view()
        assert result is not None
        assert result.n_points == 10


def test_export_2d_png_success(window, mock_file_dialog, mock_message_box):
    """Test export_2d_png success path."""
    window.state_manager.data.atoms = [MagicMock()]
    mock_file_dialog.getSaveFileName.return_value = (
        "/path/to/image.png",
        "PNG Files (*.png)",
    )
    mock_message_box.question.return_value = (
        mock_message_box.StandardButton.Yes
    )

    with (
        patch("moleditpy.ui.export_logic.AtomItem", DummyAtomItem),
        patch("moleditpy.ui.export_logic.BondItem", DummyBondItem),
        patch("moleditpy.ui.export_logic.QImage") as mock_qimage_cls,
        patch("moleditpy.ui.export_logic.QPainter") as mock_painter_cls,
    ):
        mock_image = MagicMock()
        mock_qimage_cls.return_value = mock_image
        mock_image.save.return_value = True

        mock_painter = MagicMock()
        mock_painter_cls.return_value = mock_painter
        mock_painter.begin.return_value = True
        mock_painter.isActive.return_value = True

        atom_instance = DummyAtomItem()
        window.init_manager.scene.items.return_value = [atom_instance]

        window.export_2d_png()
        mock_image.save.assert_called_with("/path/to/image.png", "PNG")
        window.statusBar().showMessage.assert_any_call(
            "2D view exported to /path/to/image.png"
        )


def test_export_2d_svg_success(window, mock_file_dialog, mock_message_box):
    """Test export_2d_svg success path."""
    window.state_manager.data.atoms = [MagicMock()]
    mock_file_dialog.getSaveFileName.return_value = (
        "/path/to/image.svg",
        "SVG Files (*.svg)",
    )
    mock_message_box.question.return_value = (
        mock_message_box.StandardButton.No
    )

    with (
        patch("moleditpy.ui.export_logic.QSvgGenerator") as mock_svg_gen_cls,
        patch("moleditpy.ui.export_logic.QPainter") as mock_painter_cls,
        patch("moleditpy.ui.export_logic.AtomItem", DummyAtomItem),
        patch("moleditpy.ui.export_logic.BondItem", DummyBondItem),
    ):
        mock_svg = MagicMock()
        mock_svg_gen_cls.return_value = mock_svg

        mock_painter = MagicMock()
        mock_painter_cls.return_value = mock_painter
        mock_painter.begin.return_value = True
        mock_painter.isActive.return_value = True

        atom_instance = DummyAtomItem()
        window.init_manager.scene.items.return_value = [atom_instance]

        window.export_2d_svg()
        mock_svg.setFileName.assert_called_with("/path/to/image.svg")
        window.statusBar().showMessage.assert_any_call(
            "2D view exported to /path/to/image.svg"
        )


def test_export_from_3d_view_no_color_logic(window):
    """Test the logic of extracting mesh without colors."""
    mock_mesh_data = TinyMesh(10)

    mock_actor = Mock(name="ActorNC")
    mock_actor.mapper = Mock()
    mock_actor.mapper.input = mock_mesh_data
    mock_actor.prop = Mock()
    mock_actor.prop.color = [1.0, 1.0, 1.0]

    window.view_3d_manager.plotter.renderer.actors = {"actor1": mock_actor}
    window.view_3d_manager.plotter.mesh = {"actor1": mock_mesh_data}

    with patch("moleditpy.ui.export_logic.pv") as mock_pv:
        mock_pv.PolyData = TinyMesh
        mock_pv.wrap.side_effect = lambda x: x

        result = window.export_from_3d_view_no_color()
        assert result is not None
        assert result.n_points == 10


def test_export_from_3d_view_with_colors_logic(window):
    """Test logic of extracting mesh with colors and splitting."""
    mock_mesh_data = TinyMesh(10)
    mock_mesh_data.point_data["colors"] = np.zeros((10, 3), dtype=np.uint8)

    mock_actor = Mock(name="ActorWC")
    mock_actor.mapper = Mock()
    mock_actor.mapper.input = mock_mesh_data
    mock_actor.prop = Mock()
    mock_actor.prop.color = [1.0, 1.0, 1.0]

    window.view_3d_manager.plotter.renderer.actors = {"actor1": mock_actor}

    with patch("moleditpy.ui.export_logic.pv") as mock_pv:
        mock_pv.PolyData = TinyMesh
        mock_pv.wrap.side_effect = lambda x: x

        results = window.export_from_3d_view_with_colors()
        assert isinstance(results, list)
        assert len(results) >= 1


def test_create_multi_material_obj_complex_cells(window):
    """Test create_multi_material_obj with Triangle Strips and Quads."""
    mock_mesh = MagicMock()
    mock_mesh.n_points = 10
    mock_mesh.points = np.zeros((10, 3))
    mock_mesh.n_cells = 2

    cell_strip = MagicMock()
    cell_strip.type = 6
    cell_strip.point_ids = [0, 1, 2, 3]

    cell_quad = MagicMock()
    cell_quad.type = 9
    cell_quad.point_ids = [4, 5, 6, 7]

    mock_mesh.get_cell.side_effect = [cell_strip, cell_quad]

    meshes = [
        {
            "mesh": mock_mesh,
            "color": [0, 255, 0],
            "name": "complex_mesh",
            "type": "actor",
            "actor_name": "complex",
        }
    ]

    with patch("builtins.open", mock_open()) as mock_file:
        window.create_multi_material_obj(meshes, "complex.obj", "complex.mtl")
        handle = mock_file()
        content = "".join([args[0] for args, _ in handle.write.call_args_list])

        assert "f 5 6 7 8" in content
        assert "f 1 2 3" in content
        assert "f 3 2 4" in content
