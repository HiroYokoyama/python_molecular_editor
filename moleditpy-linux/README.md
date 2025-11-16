# MoleditPy — A Python Molecular Editor

This is the Linux version of MoleditPy. The Open Babel fallback is disabled due to compatibility reasons.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17268532.svg)](https://doi.org/10.5281/zenodo.17268532)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

**MoleditPy** is a cross-platform, intuitive molecular editor built in Python. It provides a seamless workflow for drawing 2D molecular structures, visualizing them in 3D, and performing interactive geometric manipulations. Its powerful editing and export capabilities make it an ideal tool for preparing input files for DFT calculation software.

**Author**: HiroYokoyama
**License**: Apache-2.0
**Repository**: [https://github.com/HiroYokoyama/python\_molecular\_editor](https://github.com/HiroYokoyama/python_molecular_editor)

-----
![](img/icon.png)
![](img/screenshot.png)
-----

## Overview

This application combines a modern GUI built with **PyQt6**, powerful cheminformatics capabilities from **RDKit**, and high-performance 3D rendering by **PyVista** to provide an easy-to-use tool for drawing and visually inspecting molecular structures.

## Key Features

　Please refer to the [user manual](https://hiroyokoyama.github.io/python_molecular_editor/manual/manual) for details.

### 1\. 2D Drawing and Editing

  * **Intuitive Operations:** Easily add, edit, and delete atoms and bonds with simple mouse controls. Left-click and drag to create, and right-click to delete.
  * **Advanced Templates:** Place templates for benzene or 3- to 9-membered rings with a live preview. Features intelligent logic to automatically adjust the double-bond configuration when fusing a benzene ring to an existing structure.
  * **Charges & Radicals:** Quickly set formal charges (`+`/`-`) and radicals (`.`) on any atom using keyboard shortcuts or the toolbar.
  * **Full Element Support:** Select any element from a built-in periodic table dialog.
  * **Clipboard Operations:** Full support for `Cut` (Ctrl+X), `Copy` (Ctrl+C), and `Paste` (Ctrl+V) for molecular fragments.

### 2\. 3D Visualization and Advanced Editing

  * **2D-to-3D Conversion:** Generate a 3D conformation from your 2D drawing using RDKit's powerful algorithms. Includes a robust fallback to Open Babel if the primary method fails.
  * **Interactive 3D Editing:** Perform **comprehensive geometric editing**, including dragging atoms directly in the 3D view to intuitively create specific conformations.
  * **Precise Geometric Control:** A **powerful suite of tools** for precise geometric control:
      * **Translation:** Move the entire molecule or selected atoms to specific coordinates.
      * **Alignment:** Align the molecule by placing two selected atoms along the X, Y, or Z axis.
      * **Planarization:** Force a selection of atoms to lie on a specified plane (XY, XZ, or YZ) or project them onto their **best-fit plane**.
      * **Mirror:** Create a mirror image of the molecule across a selected plane (XY, XZ, or YZ).
      * **Set Bond Length, Angle, & Dihedral Angle:** Set precise numerical values for distances, angles, and torsional angles.
      * **Constrained Optimization:** Perform force field optimization (MMFF/UFF) while applying fixed constraints to specific distances, angles, or dihedral angles.
  * **3D Measurement Tools:** A dedicated "Measure" mode allows you to click on atoms in the 3D view to instantly measure and display:
      * **Distance** (between 2 atoms)
      * **Angle** (between 3 atoms)
      * **Dihedral Angle** (between 4 atoms)
  * **Structure Optimization:** Perform 3D structure optimization using MMFF94 or UFF force fields.
  * **Multiple Display Styles:** Visualize molecules as "Ball & Stick," "CPK (Space-filling)," "Wireframe," or "Stick" models.

### 3\. Analysis and Export

  * **Molecular Properties Analysis:** A dedicated window displays key molecular properties calculated by RDKit, such as molecular formula, weight, SMILES, LogP, and TPSA.
  * **Stereochemistry Display:** Automatically identifies and displays R/S labels for chiral centers in the 3D view after conversion.
  * **File I/O:**
      * Save and load entire sessions, including 2D/3D data and constraints, with the native `.pmeprj` project file.
      * Import structures from **MOL/SDF** files or **SMILES** strings.
      * Export 3D structures to **MOL** or **XYZ** formats, which are compatible with most DFT calculation software.
      * Export 2D and 3D views as high-resolution PNG images.

## Installation and Execution

For detailed instructions, please refer to the project [Wiki](https://github.com/HiroYokoyama/python_molecular_editor/wiki). A [Docker version](https://github.com/HiroYokoyama/python_molecular_editor_docker) is also available.

#### Requirements

`PyQt6`, `RDKit`, `NumPy`, `PyVista`, `pyvistaqt`, `openbabel`

#### Installation

1.  **Install the Package**
    This will automatically install the correct `moleditpy` (for Win/Mac) or `moleditpy-linux` (for Linux) as a dependency.

    ```bash
    pip install moleditpy-installer
    ```

2.  **Create the Shortcut**
    After installation, run this command to create the shortcut in your application menu (e.g., Start Menu or Applications folder).

    ```bash
    moleditpy-installer
    ```

-----

#### Running the Application

```bash
moleditpy
```

**(Note: The first launch may take some time while libraries like RDKit are initialized.)**

## Keyboard Shortcuts

| Key | Action | Notes |
| :--- | :--- | :--- |
| `1`/`2`/`3` | Change bond order | Single/Double/Triple bond |
| `W`/`D` | Change to stereochemical bond | Wedge / Dash bond |
| `Delete`/`Backspace` | Delete item(s) | Deletes selected or hovered items |
| `.` | Toggle radical | Cycles through 0, 1, and 2 radicals |
| `+`/`-` | Increase/Decrease charge | Changes formal charge |
| `C`, `N`, `O`, etc. | Change atom symbol | Applies to atom under cursor |
| `4` | Place benzene ring | One-shot placement on atom/bond |
| `Space` | Toggle select mode / Select all | |
| `Ctrl+J` | Perform 2D optimization (Clean Up) | |
| `Ctrl+K` | Perform 2D-to-3D conversion | |
| `Ctrl+L` | Perform 3D structure optimization | |

## Technical Details

  * **GUI and 2D Drawing (PyQt6):** The editor is built on a `QGraphicsScene`, where custom `AtomItem` and `BondItem` objects are interactively manipulated. The Undo/Redo feature is implemented by serializing the application state.
  * **Chemical Calculations (RDKit / Open Babel):** RDKit is used to generate molecule objects from 2D data, perform 3D coordinate generation, and calculate properties. Open Babel serves as a fallback for 3D conversion. All heavy computations are run on a separate `QThread` to keep the GUI responsive.
  * **3D Visualization (PyVista / pyvistaqt):** 3D rendering is achieved by generating PyVista meshes (spheres and cylinders) from RDKit conformer coordinates. A custom `vtkInteractorStyle` enables direct drag-and-drop editing of atoms in the 3D view.

## License

This project is licensed under the **Apache-2.0 License**. See the `LICENSE` file for details.

-----

# MoleditPy — A Python Molecular Editor

**MoleditPy**は、Pythonで構築されたクロスプラットフォームかつ直感的な分子エディターです。2Dでの分子描画から3D構造の可視化・編集、さらにはインタラクティブな幾何学的操作まで、シームレスなワークフローを提供します。その強力な編集機能とエクスポート機能により、**DFT計算ソフトウェアのインプット作成に最適なツール**です。

**作者**: HiroYokoyama
**ライセンス**: Apache-2.0
**リポジトリ**: [https://github.com/HiroYokoyama/python\_molecular\_editor](https://github.com/HiroYokoyama/python_molecular_editor)

-----

## 概要

このアプリケーションは、**PyQt6**によるモダンなGUI、**RDKit**による強力な化学計算、**PyVista**による高性能な3Dレンダリングを組み合わせ、分子構造の描画と視覚的な確認を容易にするツールです。

## 主な機能

詳細は、[ユーザーマニュアル](https://hiroyokoyama.github.io/python_molecular_editor/manual/manual-JP)を参照してください。

### 1\. 2D描画と編集

  * **直感的な操作:** シンプルなマウス操作で原子や結合を簡単に追加、編集、削除できます。左クリック＆ドラッグで作成し、右クリックで削除します。
  * **高度なテンプレート機能:** ベンゼン環や3〜9員環のテンプレートをライブプレビューしながら配置できます。既存の構造にベンゼン環を縮環させる際には、二重結合の配置を自動的に調整するインテリジェントなロジックを備えています。
  * **電荷とラジカル:** キーボードショートカット (`+`/`-`/`.`) やツールバーを使って、任意の原子に形式電荷やラジカルを素早く設定できます。
  * **全元素対応:** 内蔵の周期表ダイアログから任意の元素を選択できます。
  * **クリップボード操作:** 分子フラグメントの`カット` (Ctrl+X)、`コピー` (Ctrl+C)、`ペースト` (Ctrl+V) に完全対応しています。

### 2\. 3D可視化と高度な編集

  * **2D-3D変換:** 描画した2D構造から、RDKitの強力なアルゴリズムを用いて3D構造を生成します。主要な手法が失敗した場合は、Open Babelによるフォールバック機能を備えています。
  * **インタラクティブ3D編集:** 3Dビュー内の原子を直接ドラッグ操作でき、分子の形状を**本格的に編集**できます。これにより、理論計算で検討したい特定の配座を直感的に作成できます。
  * **精密な幾何学制御:** **精密な幾何学制御のための強力なツール群**を提供します。
      * **平行移動:** 分子全体または選択原子群を特定の座標へ移動します。
      * **整列:** 選択した2原子をX, Y, Z軸に沿って配置します。
      * **平面化:** 選択した3つ以上の原子を特定の平面（XY, XZ, YZ）上に配置したり、**最適フィット平面に投影**したりできます。
      * **鏡像作成:** 選択した平面 (XY, XZ, YZ) に対して分子の鏡像を作成します。
      * **結合長・角度・二面角:** 原子を選択し、目標値を入力することで、距離、角度、ねじれ角を正確に設定します。
      * **制約付き最適化:** 特定の距離、角度、二面角を固定したまま、力場計算 (MMFF/UFF) による構造最適化を実行します。
  * **3D測定ツール:** 専用の「Measure」モードで3Dビュー内の原子をクリックするだけで、以下の値を即座に測定・表示します。
      * **距離** (2原子間)
      * **角度** (3原子間)
      * **二面角** (4原子間)
  * **構造最適化:** MMFF94またはUFF力場を用いて3D構造の最適化を実行できます。
  * **多彩な表示スタイル:** 分子を「ボール＆スティック」、「CPK (空間充填)」、「ワイヤーフレーム」、「スティック」モデルで表示できます。

### 3\. 解析とエクスポート

  * **分子特性解析:** 専用ウィンドウに、分子式、分子量、SMILES、LogP、TPSAなど、RDKitによって計算された主要な分子特性を一覧表示します。
  * **立体化学表示:** 3D変換後、キラル中心を自動的に認識し、R/Sラベルを3Dビューに表示します。
  * **ファイル入出力:**
      * 2D/3Dデータや制約情報を含むセッション全体を、独自のプロジェクトファイル (`.pmeprj`) として保存・読み込みできます。
      * **MOL/SDF**ファイルや**SMILES**文字列から構造をインポートできます。
      * 3D構造を**MOL**または**XYZ**形式でエクスポートでき、これらは多くのDFT計算ソフトウェアと互換性があります。
      * 2Dおよび3Dビューを高解像度のPNG画像としてエクスポートできます。

## インストールと実行

詳細な手順については、プロジェクトの[Wiki](https://github.com/HiroYokoyama/python_molecular_editor/wiki)を参照してください。[Docker版](https://github.com/HiroYokoyama/python_molecular_editor_docker)も利用可能です。

#### 必要ライブラリ

`PyQt6`, `RDKit`, `NumPy`, `PyVista`, `pyvistaqt`, `openbabel`

#### インストール

1.  **パッケージのインストール**
    このコマンドを実行すると、お使いのOS（Windows/macOSまたはLinux）に適した `moleditpy` 本体が自動的にインストールされます。

    ```bash
    pip install moleditpy-installer
    ```

2.  **ショートカットの作成**
    インストール後、このコマンドを実行すると、アプリケーションメニュー（スタートメニューやアプリケーションフォルダなど）にショートカットが作成されます。

    ```bash
    moleditpy-installer
    ```

#### アプリケーションの起動

```bash
moleditpy
```

**（注：初回起動時は、RDKitなどのライブラリの初期化のため、起動に時間がかかる場合があります。）**

## キーボードショートカット

| キー | 操作 | 補足 |
| :--- | :--- | :--- |
| `1`/`2`/`3` | 結合次数を変更 | 単結合/二重結合/三重結合 |
| `W`/`D` | 立体化学結合に変更 | Wedge / Dash 結合 |
| `Delete` / `Backspace` | アイテムの削除 | 選択またはカーソル下のアイテムを削除 |
| `.` | ラジカルをトグル | 0, 1, 2ラジカルを循環 |
| `+`/`-` | 電荷を増減 | 形式電荷の変更 |
| `C`, `N`, `O` など | 原子記号を変更 | カーソル下の原子に適用 |
| `4` | ベンゼン環の配置 | カーソル下の原子/結合にワンショットで配置 |
| `Space` | 選択モード切替 / 全選択 | |
| `Ctrl+J` | 2D最適化を実行 | |
| `Ctrl+K` | 3D変換を実行 | |
| `Ctrl+L` | 3D最適化を実行 | |

## 技術的な仕組み

  * **GUIと2D描画 (PyQt6):** `QGraphicsScene`上にカスタムの`AtomItem`（原子）と`BondItem`（結合）を配置し、対話的に操作します。Undo/Redo機能は、アプリケーションの状態をシリアライズしてスタックに保存することで実現しています。
  * **化学計算 (RDKit / Open Babel):** 2DデータからRDKit分子オブジェクトを生成し、3D座標生成や分子特性計算を実行します。RDKitでの3D座標生成が失敗した際は、Open Babelにフォールバックします。重い計算処理は別スレッド (`QThread`) で実行し、GUIの応答性を維持しています。
  * **3D可視化 (PyVista / pyvistaqt):** RDKitのコンフォーマ座標からPyVistaのメッシュ（球や円柱）を生成して描画します。カスタムの`vtkInteractorStyle`を実装し、3Dビュー内での原子の直接的なドラッグ＆ドロップ編集を可能にしています。

## ライセンス

このプロジェクトは **Apache-2.0 License** のもとで公開されています。詳細は `LICENSE` ファイルを参照してください。
