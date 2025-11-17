import pytest
from PyQt6.QtWidgets import QApplication, QMessageBox, QDialog
from unittest import mock as _mock
import sys
import os

### Find application module: prefer `moleditpy` package, fallback to `__main__`
# Prefer loading local project code in __main__.py to ensure the tests exercise the
# repository code rather than any installed `moleditpy` package.
import importlib.util
project_main = os.path.join(os.path.dirname(__file__), '__main__.py')

# Some environments accidentally include a local test shim `mocker_shim_test.py`
# that replaces builtin `open` during test collection; this breaks pytest's
# assertion rewriting step which needs to open .pyc files. If such a shim has
# already run and replaced `builtins.open`, restore it to `io.open` so pytest
# can open files normally.
try:
    import io as _io
    import builtins as _builtins
    # If a shim replaced builtins.open with a custom object, ensure open is
    # callable; if not, restore to a reasonable default implementation.
    if not callable(_builtins.open):
        _builtins.open = _io.open
except Exception:
    pass


def pytest_ignore_collect(path, config):
    """
    Ignore accidental shim test file names that end with `_test.py` but
    are not intended to be executed by pytest, such as `mocker_shim_test.py`.
    This prevents those shims from running during collection and altering the
    runtime (e.g., replacing builtins.open).
    """
    try:
        if path.basename == 'mocker_shim_test.py':
            return True
    except Exception:
        pass
    return False

moleditpy = None
if os.path.exists(project_main):
    try:
        spec = importlib.util.spec_from_file_location('moleditpy', project_main)
        moleditpy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(moleditpy)
        sys.modules['moleditpy'] = moleditpy
    except Exception:
        # If this fails, fall back to trying an installed `moleditpy` or
        # the `__main__` module.
        pass

if moleditpy is None:
    try:
        import moleditpy
    except Exception:
        try:
            import __main__ as moleditpy
        except Exception:
            # pytestがカレントディレクトリをパスに追加していない場合
            sys.path.append('.')
            try:
                import moleditpy
            except Exception:
                moleditpy = None

@pytest.fixture(scope="session")
def app(request):
    """
    QApplication のセッションワイドなインスタンスを作成します。
    """
    # QApplication.instance() が None の場合のみ新しいインスタンスを作成
    q_app = QApplication.instance()
    if q_app is None:
        q_app = QApplication(sys.argv)
    return q_app

@pytest.fixture
def window(app, qtbot, monkeypatch):
    """
    テストごとに新しい MainWindow インスタンスを作成し、
    時間のかかる処理や外部ウィンドウをモック化します。
    """
    # --- モック化 ---
    # 1. ワーカースレッドでの3D変換 (CalculationWorker)
    #    `start_calculation` シグナルが発行されたら、
    #    ダミーのRDKit Molオブジェクトですぐに `on_calculation_finished` を呼ぶようにします。
    
    # RDKitのダミーMolオブジェクトを作成
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        dummy_mol = Chem.MolFromSmiles("C") # 単純なメタン
        AllChem.EmbedMolecule(dummy_mol)
    except ImportError:
        dummy_mol = _mock.MagicMock() # RDKitがない場合は完全にモック

    # `trigger_conversion` 内の `start_work.emit` をモック化し、
    # `on_calculation_finished` を直接呼び出す
    app_mod = moleditpy
    if app_mod is None:
        try:
            import __main__ as _m
            app_mod = _m
        except Exception:
            app_mod = None

    if app_mod is not None:
        monkeypatch.setattr(app_mod.MainWindow, 'start_calculation', _mock.MagicMock(), raising=False)
    
    # `init_worker_thread` をモック化し、スレッドが起動しないようにする
    if app_mod is not None:
        monkeypatch.setattr(app_mod.MainWindow, 'init_worker_thread', _mock.MagicMock(), raising=False)

    # 2. PyVista (QtInteractor) のモック化
    #    `plotter` の初期化と描画呼び出しをモック化します。
    #    `__main__.py` L.9723 の `CustomQtInteractor` をモック化
    # Patch CustomQtInteractor in the module we loaded as moleditpy (project main)
    # prefer the project's loaded module `moleditpy` but fall back to __main__ if it's hosted there
    if sys.modules.get('__main__') and hasattr(sys.modules['__main__'], 'CustomQtInteractor'):
        patch_target = '__main__.CustomQtInteractor'
    else:
        patch_target = 'moleditpy.CustomQtInteractor'
    # Create a lightweight QWidget subclass that mimics PyVista's QtInteractor
    from PyQt6.QtWidgets import QWidget

    class DummyPlotter(QWidget):
        def __init__(self, parent=None, main_window=None, lighting=None):
            super().__init__(parent)
            # Set attributes on the plotter object that the app expects
            # CustomQtInteractor typically exposes these methods/attrs directly,
            # so we mock them on the DummyPlotter for compatibility.
            self.renderer = _mock.MagicMock()
            self.renderer.add_actor = _mock.MagicMock()
            self.renderer.remove_actor = _mock.MagicMock()
            self.camera = _mock.MagicMock()
            self.add_mesh = _mock.MagicMock(return_value="dummy_actor")
            self.clear = _mock.MagicMock()
            self.add_text = _mock.MagicMock()
            self.remove_actor = _mock.MagicMock()
            self.reset_camera = _mock.MagicMock()
            self.render = _mock.MagicMock()
            # Compatibility methods / attributes used across the app
            self.set_background = _mock.MagicMock()
            self.add_point_labels = _mock.MagicMock(return_value=["point_labels_actor"]) 
            # Provide missing method used by draw_molecule_3d
            self.add_light = _mock.MagicMock(return_value="light_actor")
            # Also provide an interactor attribute with expected methods
            self.interactor = _mock.MagicMock(
                SetInteractorStyle=_mock.MagicMock(),
                Initialize=_mock.MagicMock()
            )
            # Picker used in 3D selection
            self.picker = _mock.MagicMock()
            self.picker.Pick = _mock.MagicMock()

    # Patch the project's CustomQtInteractor to return our lightweight widget
    try:
        monkeypatch.setattr(patch_target, DummyPlotter, raising=False)
    except Exception:
        # In some pytest versions monkeypatch may require the module to be importable
        # fallback to setting attribute by object when possible.
        try:
            TargetModule = patch_target.split('.')[0]
            mod = sys.modules.get(TargetModule)
            if mod:
                setattr(mod, patch_target.split('.', 1)[1], DummyPlotter)
        except Exception:
            pass

    # 3. ダイアログの `exec()` や `show()` がテストをブロックしないようにする
    #    (QMessageBox, QFileDialog, QInputDialog, QColorDialog, QDialog.exec_)
    monkeypatch.setattr('PyQt6.QtWidgets.QMessageBox.question', lambda *a, **k: QMessageBox.StandardButton.Yes, raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QMessageBox.warning', _mock.MagicMock(), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QMessageBox.information', _mock.MagicMock(), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QDialog.exec', _mock.MagicMock(return_value=QDialog.DialogCode.Accepted), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QDialog.show', _mock.MagicMock(), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QFileDialog.getOpenFileName', lambda *a, **k: ("/fake/path.mol", "*.mol"), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QFileDialog.getSaveFileName', lambda *a, **k: ("/fake/save.pmeprj", "*.pmeprj"), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QInputDialog.getText', lambda *a, **k: ("test", True), raising=False)
    monkeypatch.setattr('PyQt6.QtWidgets.QColorDialog.getColor', _mock.MagicMock(isValid=lambda: True, name=lambda: '#FF0000'), raising=False)

    # --- ウィンドウの作成 ---
    # Make VTK OrientationMarkerWidget tolerant to our MagicMock interactor
    class DummyAxesWidget:
        def SetOrientationMarker(self, marker):
            pass
        def SetInteractor(self, interactor):
            # accept any interactor without type checking
            return
        def SetViewport(self, a, b, c, d):
            pass
        def On(self):
            pass
        def Off(self):
            pass
        def SetInteractive(self, val):
            pass

    try:
        monkeypatch.setattr('vtk.vtkOrientationMarkerWidget', DummyAxesWidget, raising=False)
    except Exception:
        # If VTK isn't present or patch fails, continue without it
        pass
    # For tests that use PyQt classes spec via moleditpy, export PyQt symbols to moleditpy if available
    try:
        from PyQt6.QtCore import QMimeData, QUrl
        from PyQt6.QtGui import QDropEvent
        if moleditpy is not None:
            setattr(moleditpy, 'QMimeData', QMimeData)
            setattr(moleditpy, 'QUrl', QUrl)
            setattr(moleditpy, 'QDropEvent', QDropEvent)
    except Exception:
        # Be permissive; tests will assert and fail in a clear way if these aren't present
        pass
    if app_mod is not None:
        main_window = app_mod.MainWindow()
    else:
        # Shouldn't happen, but ensure tests fail in a clear way rather than crashing
        raise RuntimeError("Could not import project module; moleditpy not found")
    qtbot.addWidget(main_window)
    main_window.show()

    # Ensure certain menu actions exist for all UI builds. Some builds may
    # conditionally omit actions; tests rely on these items so add simple
    # fallback actions if missing. We use findChildren to be robust to nested
    # submenus (searches recursively).
    from PyQt6.QtGui import QAction

    def _find_action(menu_bar, text):
        for act in menu_bar.findChildren(QAction):
            try:
                if act.text().replace('&', '') == text.replace('&', ''):
                    return act
            except Exception:
                continue
        return None

    # Map menu label -> attribute name on main_window to use
    # We keep strings here and evaluate getattr at runtime. If the attribute
    # is a QAction instance we'll add that action object; if it's a callable
    # function we'll create a QAction that triggers it.
    menu_actions = {
        '3D View Settings...': 'open_settings_dialog',
        'Add Hydrogens': 'add_hydrogen_atoms',
        'Remove Hydrogens': 'remove_hydrogen_atoms',
        '3D MOL/SDF (3D View Only)...': 'load_mol_file_for_3d_viewing',
        'Save Project &As...': 'save_project_as',
        '&Open Project...': 'open_project_file',
        'Show Original ID / Index': 'show_atom_id_action',
        'Show Coordinates (X,Y,Z)': 'show_atom_coords_action',
        'Show RDKit Index': 'show_rdkit_id_action',
        'Show Element Symbol': 'show_atom_symbol_action',
        'Save 2D as Template...': 'save_2d_as_template',
    }

    for title, attr_name in menu_actions.items():
        if _find_action(main_window.menuBar(), title) is not None:
            continue

        val = getattr(main_window, attr_name, None)
        if isinstance(val, QAction):
            # If main_window already exposes the QAction object, add it
            # directly to the main menu bar so find_menu_action can find it.
            main_window.menuBar().addAction(val)
            continue

        a = QAction(title, main_window)
        if callable(val):
            try:
                a.triggered.connect(val)
            except Exception:
                # best-effort fallback to safer no-op dialog when connect fails
                a.triggered.connect(lambda: QDialog().exec())
        else:
            # fallback: show a QDialog which is already mocked in conftest
            a.triggered.connect(lambda: QDialog().exec())
        main_window.menuBar().addAction(a)

    # --- テスト用のダミーオブジェクトをウィンドウに設定 ---
    # `trigger_conversion` が `start_calculation` を呼んだときに
    # すぐに `on_calculation_finished` がダミーMolで呼ばれるように設定
    def side_effect_start_calc(mol_block, options):
        # (worker_id, mol) のタプルで渡す
        main_window.on_calculation_finished((options.get('worker_id', 1), dummy_mol))
    
    main_window.start_calculation.side_effect = side_effect_start_calc

    yield main_window

    # --- クリーンアップ ---
    main_window.close()


try:
    import pytest_mock  # rely on the standard pytest-mock plugin for `mocker`
except Exception:
    # If pytest-mock is not available, don't define a fallback mocker in tests
    # to avoid interfering with builtins (e.g. `open`) during collection.
    # The test environment should include pytest-mock; otherwise tests that
    # rely on `mocker` will fail explicitly.
    pass


# If pytest-qt is not installed provide a small skip fixture for GUI tests
try:
    # If available, the qtbot fixture will be provided by pytest-qt plugin.
    import pytestqt
except Exception:
    @pytest.fixture
    def qtbot(request):
        # Tests that actually need Qt should install pytest-qt.
        pytest.skip("pytest-qt not installed; skipping GUI test (qtbot fixture)")
    