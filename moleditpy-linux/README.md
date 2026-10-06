# MoleditPy — A Python Molecular Editor

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17268532.svg)](https://doi.org/10.5281/zenodo.17268532)
[![PyPI version](https://badge.fury.io/py/MoleditPy.svg)](https://badge.fury.io/py/MoleditPy)
[![Python Versions](https://img.shields.io/badge/python-3.9--3.14-blue.svg)](https://pypi.org/project/MoleditPy/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://github.com/HiroYokoyama/python_molecular_editor/actions/workflows/tests.yml/badge.svg)](https://github.com/HiroYokoyama/python_molecular_editor/actions)
[![codecov](https://codecov.io/gh/HiroYokoyama/python_molecular_editor/graph/badge.svg)](https://codecov.io/gh/HiroYokoyama/python_molecular_editor)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/moleditpy?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/moleditpy)
[![](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/HiroYokoyama)

**MoleditPy** is a programmable, cross-platform molecular editor built with PyQt6, RDKit, and PyVista. It takes you from a 2D sketch to an editable 3D structure and out to input files for DFT calculations, and it can be extended with plain Python plugins.

**[Website](https://hiroyokoyama.github.io/python_molecular_editor/)** · **[User manual](https://hiroyokoyama.github.io/python_molecular_editor/manual/manual)** ([日本語](https://hiroyokoyama.github.io/python_molecular_editor/manual/manual-JP)) · **[Plugins](https://hiroyokoyama.github.io/moleditpy-plugins/)** · **[Wiki](https://github.com/HiroYokoyama/python_molecular_editor/wiki)**

![MoleditPy screenshot](https://hiroyokoyama.github.io/python_molecular_editor/img/screenshot.png)

## Features

- **2D drawing** — atoms, bonds, ring templates, alkyl chains, charges, radicals, stereo bonds, and the full periodic table.
- **3D editing** — 2D-to-3D conversion, dragging atoms in the 3D view, exact bond lengths, angles, and dihedrals, and MMFF94/UFF optimization with constraints.
- **Analysis and export** — molecular properties, R/S labels, MOL/SDF/XYZ/SMILES/InChI import and export, PNG/SVG images, and STL/OBJ for 3D printing.
- **Plugins** — drop a Python file into the plugin folder to add menu actions and tools. Browse the official collection in the [Plugin Explorer](https://hiroyokoyama.github.io/moleditpy-plugins/explorer/).

See the [user manual](https://hiroyokoyama.github.io/python_molecular_editor/manual/manual) for every feature and keyboard shortcut.

## Installation

```bash
pip install moleditpy-installer
python -m moleditpy_installer
```

The first command installs the right package for your platform (`moleditpy`, or `moleditpy-linux` on Linux); the second creates an application shortcut. Start the app from the shortcut or with:

```bash
moleditpy
```

The first launch can take a while as RDKit and the other libraries initialize. A [Windows installer](https://hiroyokoyama.github.io/python_molecular_editor/windows-installer/windows_installer), a [macOS app bundle](https://hiroyokoyama.github.io/python_molecular_editor/macos-installer/macos_installer), and a [Docker image](https://github.com/HiroYokoyama/python_molecular_editor_docker) are also available.

> **Security note:** the legacy `.pmeraw` project format uses Python pickle, which can execute arbitrary code when loaded. Only open `.pmeraw` files you created yourself, and use the JSON-based `.pmeprj` format for sharing.

## Citation

If you use this software in your work, please cite it as follows:

```
Yokoyama, H. (2026). MoleditPy — A Python-based molecular editing software. Zenodo. https://doi.org/10.5281/zenodo.17268532
```

Please also cite the plugins you used; see [Citation in the plugin collection](https://github.com/HiroYokoyama/moleditpy-plugins#citation).

## License & Disclaimer

This project is licensed under the GNU General Public License v3.0 (GPLv3) - see the [LICENSE](https://github.com/HiroYokoyama/python_molecular_editor/blob/main/LICENSE) file for details. As open-source software, it is provided 'as is' without warranty of any kind, and the author assumes no responsibility or liability for the results. Although outputs have been carefully verified, users are strongly encouraged to independently check and validate them for critical applications (such as publications). If you encounter any bugs, please open an issue.
