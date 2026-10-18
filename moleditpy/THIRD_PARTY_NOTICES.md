# Third-party notices

MoleditPy is LGPL-3.0-or-later. Its license does not replace third-party licenses.

| Component | License / source |
| --- | --- |
| PySide6-Essentials and Shiboken6 | LGPL-3.0 (or other terms offered by Qt); https://code.qt.io/pyside/pyside-setup.git |
| Qt Core, GUI, Widgets, SVG and required support libraries | LGPL-3.0 (or other terms offered by Qt); https://code.qt.io/qt/qtbase.git and https://code.qt.io/qt/qtsvg.git |
| RDKit | BSD-3-Clause; https://github.com/rdkit/rdkit |
| PyVista and pyvistaqt | MIT; https://github.com/pyvista/pyvista and https://github.com/pyvista/pyvistaqt |
| VTK | BSD-3-Clause; https://gitlab.kitware.com/vtk/vtk |
| NumPy | BSD-3-Clause and bundled notices; https://github.com/numpy/numpy |
| Python | PSF license and included notices; https://github.com/python/cpython |
| PyInstaller bootloader | GPL-2.0-or-later with the PyInstaller exception; https://github.com/pyinstaller/pyinstaller |

Dependency wheels and frozen distributions include additional libraries (for example image codecs, numerical libraries, and platform support). Preserve their accompanying license and copyright notices. Build archives include installed dependency license directories under `licenses/third-party` and a version/source manifest.

The core does not depend on or bundle Open Babel. Optional plugins and their dependencies carry independent license notices. The frozen build deliberately excludes unused Qt modules, including GPL-only add-on modules; installing arbitrary Qt add-ons can change a distributor's obligations.

Official Qt licensing information: https://doc.qt.io/qt-6/licensing.html and https://doc.qt.io/qtforpython-6/licenses.html. GNU license texts: https://www.gnu.org/licenses/lgpl-3.0.html and https://www.gnu.org/licenses/gpl-3.0.html.
