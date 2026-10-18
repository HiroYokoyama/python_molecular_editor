# Distribution, source, and replacing libraries

MoleditPy is LGPL-3.0-or-later. Read `LICENSE` together with `COPYING`. Recipients may inspect and modify the program, replace LGPL libraries, and reverse engineer for debugging their modifications. No warranty is provided.

## Source distributions and installed Python packages

The corresponding MoleditPy source, tests, build specifications and release workflows are in https://github.com/HiroYokoyama/python_molecular_editor. Each frozen build includes `licenses/build-manifest.json` with its exact source commit and dependency versions. Release archives must accompany or link the exact corresponding source and build instructions. The wheel and source archive include both GNU license texts and these notices.

For a Python installation, create a virtual environment, install `moleditpy` and substitute compatible PySide6-Essentials/Shiboken6 wheels in that environment with pip. Running from source uses those installed shared libraries directly. Example:

```sh
python -m venv .venv
# Activate .venv with the command for your shell.
python -m pip install -e ./moleditpy
python -m pip install --force-reinstall /path/to/modified/shiboken6.whl /path/to/modified/PySide6_Essentials.whl
python -m moleditpy --safe
```

## Frozen Windows and macOS applications

The application is an unpacked PyInstaller directory / `.app` bundle, with Qt and PySide6 shared libraries stored outside the executable archive. It does not statically link those libraries into MoleditPy. Keep a backup, close the application, and replace a compatible set of Qt libraries and matching PySide6/Shiboken6 bindings. Use the architecture, Qt ABI and versions recorded in `licenses/build-manifest.json`.

On Windows the files are normally below `MoleditPy/_internal/PySide6/` (Qt6Core.dll, Qt6Gui.dll, Qt6Widgets.dll and associated bindings), with Shiboken6 below `_internal/shiboken6/`. Replace the corresponding files and plugin directories while preserving the layout. On macOS inspect `MoleditPy.app/Contents/Frameworks/PySide6/` and `Contents/Frameworks/`; replace matching `.framework`, `.dylib` and binding files and retain relative paths. Bundled paths vary with PyInstaller versions; locate the external `QtCore`/`Qt6Core` shared library first. macOS may require ad-hoc signing of a locally modified bundle: `codesign --force --deep --sign - MoleditPy.app`. This is a local platform requirement; MoleditPy does not require the publisher's signature.

A reliable alternative to file replacement is to rebuild with modified dependencies: check out the manifest commit, install the modified compatible wheels in a clean environment, install `moleditpy` and PyInstaller, then run `pyinstaller --clean --noconfirm windows-installer/script/MoleditPy.spec` or `macos-installer/script/MoleditPy.spec` from the repository root. The build uses normal dynamic Qt libraries and does not impose a signature or activation lock.

## Corresponding dependency source

The manifest records installed versions and upstream source locations. Obtain Qt/PySide6 sources matching those versions from https://download.qt.io/official_releases/QtForPython/ and https://download.qt.io/official_releases/qt/, or their matching tags at https://code.qt.io/pyside/pyside-setup.git and https://code.qt.io/qt/qt5.git. Follow that version's upstream build instructions and preserve its third-party source notices. Distribution of modified libraries also requires their modified corresponding source under their applicable licenses. Release distributors must make required corresponding source available alongside the binaries, or use another compliant LGPL/GPL source-delivery method; an upstream project link alone is not a substitute for that obligation.
