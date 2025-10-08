# MoleditPy -- Python Molecular Editor

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17268532.svg)](https://doi.org/10.5281/zenodo.17268532)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

A cross-platform, simple, and intuitive molecular structure editor built in Python. It allows 2D molecular drawing and 3D structure visualization. It supports exporting structure files, and its interactive 3D editing capabilities allow for the intuitive creation of desired conformations, making it **ideal for preparing inputs for DFT calculation software**.

**Author**: HiroYokoyama
**License**: Apache-2.0
**Repository**: [https://github.com/HiroYokoyama/python\_molecular\_editor](https://github.com/HiroYokoyama/python_molecular_editor)

-----

Pythonで構築された、クロスプラットフォームかつシンプルで直感的な分子構造エディターです。2Dでの分子描画と3D構造可視化ができます。構造ファイルのエクスポートをサポートし、さらにインタラクティブな3D編集機能で目的の配座を直感的に作成できるため、**DFT計算のインプット作成に最適です**。

**作者**: HiroYokoyama
**ライセンス**: Apache-2.0
**リポジトリ**: [https://github.com/HiroYokoyama/python\_molecular\_editor](https://github.com/HiroYokoyama/python_molecular_editor)

-----
![](img/screenshot.png)
-----

## Overview

This application is a tool for easily drawing molecular structures and visually inspecting their three-dimensional shapes. It combines a modern GUI by **PyQt6**, powerful chemical calculations by **RDKit**, and high-performance 3D rendering by **PyVista**.

-----

### 概要

このアプリケーションは、分子構造を容易に描き、その3次元的な形状を視覚的に確認するためのツールです。**PyQt6**によるモダンなGUI、**RDKit**による強力な化学計算、**PyVista**による高性能な3Dレンダリングを組み合わせています。

-----

## Key Features

### 1\. 2D Drawing and Editing

  * **Intuitive Operations:** Easily add, edit, or delete atoms and bonds with intuitive mouse controls. Add/edit with click-and-drag, and **delete items with a right-click** or by using the `Delete` / `Backspace` keys.
  * **Improved Template Placement:** Place templates for benzene or 3- to 9-membered rings with a live preview. Features **advanced logic to automatically adjust the double-bond configuration of benzene rings** when fused to existing bonds.
  * **Charge/Radical Operations:** Easily set formal charges and radicals by clicking on an atom or using keyboard shortcuts (`+`/`-`/`.`).
  * **Selection from Periodic Table:** Select any element from a periodic table dialog.
  * **Enhanced Selection Tools:**
      * Supports `Cut` (`Ctrl+X`), `Copy` (`Ctrl+C`), and `Paste` (`Ctrl+V`) for clipboard operations on molecular fragments.
      * `Space`: Toggles select mode / Selects all in select mode.


### 1\. 2D描画と編集

  * **直感的な操作:** 直感的なマウス操作で原子や結合の追加・編集・削除が簡単に行えます。クリック＆ドラッグで追加・編集し、**右クリックでアイテムを削除**できます。また`Delete` / `Backspace`キーでも削除可能です。
  * **テンプレート配置の改良:** ベンゼン環や3〜9員環のテンプレートをプレビューして配置可能。既存の結合にフューズする際、**ベンゼン環の二重結合配置を自動調整する高度なロジック**を搭載しました。
  * **電荷・ラジカル操作:** 原子をクリックするか、キーボードショートカット (`+`/`-`/`.`) で、形式電荷やラジカルを簡単に設定できます。
  * **周期表からの選択:** 周期表ダイアログから任意の元素を選択可能です。
  * **選択操作の充実:**
      * `Cut` (`Ctrl+X`), `Copy` (`Ctrl+C`), `Paste` (`Ctrl+V`) に対応し、分子フラグメントのクリップボード操作が可能です。
      * `Space`: 選択モード切替 / 選択モードで全選択。

-----

### 2\. Keyboard Shortcuts

| Key | Action | Notes |
| :--- | :--- | :--- |
| `1`/`2`/`3` | Change bond order | Single/Double/Triple bond |
| `W`/`D` | Change to stereochemical bond | Wedge / Dash bond |
| `Delete` / `Backspace` | Delete item(s) | Deletes selected or hovered items |
| `.` | Toggle radical | 0 → 1 → 2 → 0 |
| `+`/`-` | Increase/Decrease charge | Change formal charge |
| `C`, `N`, `O`, etc. | Change atom symbol | Applies to atom under cursor or selection |
| `4` | Place benzene ring | One-shot placement on atom/bond under cursor |
| `Ctrl+J` | Perform 2D optimization | |
| `Ctrl+K` | Perform 3D conversion | |
| `Ctrl+L` | Perform 3D optimization | |


### 2\. キーボードショートカット

| キー | 操作 | 補足 |
| :--- | :--- | :--- |
| `1`/`2`/`3` | 結合次数を変更 | 単結合/二重結合/三重結合 |
| `W`/`D` | 立体化学結合に変更 | Wedge / Dash 結合 |
| `Delete` / `Backspace` | アイテムの削除 | 選択またはカーソル下のアイテムを削除 |
| `.` | ラジカルをトグル | 0 → 1 → 2 → 0 |
| `+`/`-` | 電荷を増減 | 形式電荷の変更 |
| `C`, `N`, `O` など | 原子記号を変更 | カーソル下または選択中の原子に適用 |
| `4` | ベンゼン環の配置 | カーソル下の原子/結合にワンショットで配置 |
| `Ctrl+J` | 2D最適化を実行 | |
| `Ctrl+K` | 3D変換を実行 | |
| `Ctrl+L` | 3D最適化を実行 | |

-----

### 3\. 2D Structure Optimization

  * Performs automatic layout using RDKit's `Compute2DCoords` (**Optimize 2D**).
  * Implements logic to detect and **automatically separate and resolve overlapping groups of atoms** (such as non-bonded fragments).


### 3\. 2D構造の最適化

  * RDKit の `Compute2DCoords` を使った自動レイアウト（**Optimize 2D**）を実行します。
  * 完全に重なった原子グループ（結合していないフラグメント同士など）を検出し、**自動で分離・解消するロジック**を実装しました。

-----

### 4\. High-Quality 3D Visualization and Analysis

  * **3D Conversion:** Generates 3D coordinates with RDKit and optimizes them using an MMFF94-based force field (**Convert to 3D**). If RDKit fails, it executes a **fallback to Open Babel** to enhance robustness.
  * **Interactive Display:** Provides interactive 3D visualization (Ball & Stick / CPK styles) powered by PyVista / pyvistaqt.
  * **Interactive 3D Editing:** When 3D edit mode is enabled, you can **directly drag atoms in the 3D view with the mouse** to fine-tune their positions. This allows for the intuitive creation of specific conformations you wish to investigate in theoretical calculations.
  * **Chiral Label Display:** After 3D conversion, automatically assigns R/S labels to chiral centers and displays them in the **3D view**.
  * **Molecule Analysis Window:** A dedicated window that **lists key molecular properties** based on RDKit, such as molecular formula, molecular weight, SMILES, LogP, and TPSA.


### 4\. 高品質な3D可視化と分析

  * **3D変換:** RDKit で 3D 座標を生成し MMFF94 ベースで最適化（**Convert to 3D**）します。RDKitでの生成に失敗した場合、**Open Babelによるフォールバック**を実行し、堅牢性を高めています。
  * **インタラクティブ表示:** PyVista / pyvistaqt によるインタラクティブな3D表示（Ball & Stick / CPK スタイル）を提供します。
  * **インタラクティブ3D編集:** 3D編集モードを有効にすると、3Dビュー内の原子を**マウスで直接ドラッグして位置を微調整**できます。これにより、理論計算で検討したい特定の配座（コンフォメーション）を直感的に作り出すことができます。
  * **キラルラベル表示:** 3D変換後、キラル中心に R/S ラベルを自動で付与し、**3Dビュー**に表示します。
  * **分子分析ウィンドウ:** 分子式、分子量、SMILES、LogP、TPSAなど、RDKitベースの**主要な分子特性を一覧表示**する専用ウィンドウがあります。

-----

### 5\. File I/O

  * **Project File (.pmeraw):** Allows you to completely save and load the editing state, including 2D drawing data, 3D structure, charges, radicals, and chiral label states.
  * Save 2D structures in MOL format.
  * Save 3D structures in MOL / XYZ formats. These formats are widely supported as inputs for many DFT calculation software packages.
  * Supports importing from MOL/SDF files.
  * Supports importing from a SMILES string.


### 5\. ファイル入出力

  * **プロジェクトファイル (.pmeraw):** 2D描画データと3D構造、電荷、ラジカル、キラルラベルの状態など、編集状態を完全に保存/読み込みできます。
  * 2D を MOL 形式で保存。
  * 3D を MOL / XYZ 形式で保存。これらの形式は多くのDFT計算ソフトウェアでインプットとして利用できます。
  * MOL/SDF の読み込みに対応しています。
  * SMILES文字列からのインポートに対応しています。

-----

## Execution and Installation

#### Requirements

`PyQt6`, `RDKit`, `NumPy`, `PyVista`, `pyvistaqt`, `openbabel`

#### Installation Example

**Using pip:**

```bash
pip install moleditpy
```

> **Note**
> It is recommended to install RDKit and Open Babel using `conda`, a scientific computing distribution. Open Babel is required for the 3D structure conversion fallback feature.

#### Running the App

```bash
moleditpy
```


## 実行とインストール

#### 必要ライブラリ

`PyQt6`, `RDKit`, `NumPy`, `PyVista`, `pyvistaqt`, `openbabel`

#### インストール例

**pip を使う場合:**

```bash
pip install moleditpy
```

> **Note**
> RDKit と Open Babel は、科学計算ディストリビューションである `conda` を使ってインストールすることが推奨されます。Open Babelは3D構造変換のフォールバック機能に必要です。

#### アプリの起動

```bash
moleditpy
```

-----

## Development Environment

The recommended development and execution environment for MoleditPy is as follows.

| Component | Version |
| :--- | :--- |
| **Python** | `3.13.7` |
| **numpy** | `2.3.3` |
| **openbabel-wheel** | `3.1.1.22` |
| **PyQt6** | `6.9.1` |
| **QtPy** | `2.4.3` |
| **rdkit** | `2025.3.6` |
| **vtk** | `9.5.2` |
| **pyvista** | `0.46.3` |
| **pyvistaqt** | `0.11.3` |


### 開発環境

MoleditPyの動作確認および推奨される開発・実行環境は以下の通りです。

| コンポーネント | バージョン |
| :--- | :--- |
| **Python** | `3.13.7` |
| **numpy** | `2.3.3` |
| **openbabel-wheel** | `3.1.1.22` |
| **PyQt6** | `6.9.1` |
| **QtPy** | `2.4.3` |
| **rdkit** | `2025.3.6` |
| **vtk** | `9.5.2` |
| **pyvista** | `0.46.3` |
| **pyvistaqt** | `0.11.3` |

-----

## Windows: Shortcut and File Association Guide

This guide provides detailed steps for setting up the application in a Windows environment, allowing you to open project files (`.pmeraw`) by double-clicking and launch the application from a shortcut with a custom icon.

### 1\. Paths for Executable and Icon

| Item | Details |
| :--- | :--- |
| **Target Application** | `C:\Users\%USERNAME%\AppData\Local\Programs\Python\PythonXX\Scripts\moleditpy.exe` |
| **File Extension** | `.pmeraw` |
| **Icon File Path** | `C:\Users\%USERNAME%\AppData\Local\Programs\Python\PythonXX\Lib\site-packages\moleditpy\assets\icon.ico` |

> **Note:** Replace `%USERNAME%` with your actual username and `PythonXX` with the directory name of your installed Python version (e.g., `Python311`).

### 2\. File Association Steps (to open .pmeraw files on double-click)

Follow these steps to associate `.pmeraw` files with `moleditpy.exe`.

1.  **Find a `.pmeraw` file:** In File Explorer, **right-click** on a file you want to associate (e.g., `sample.pmeraw`).
2.  **Change association:** From the context menu, select **"Properties"**, and in the "General" tab, click the **"Change..."** button next to "Opens with:".
3.  **Select the program:** In the "How do you want to open this file?" window, choose "More apps" → **"Look for another app on this PC"**.
4.  **Specify the executable:** Navigate to the path of **`moleditpy.exe`** (with placeholders replaced) mentioned above and click "Open" to complete the process.

> **Tip:** Applying a custom icon to a file association system-wide typically requires editing the registry. This is an advanced operation and should be done with caution.

### 3\. Desktop Shortcut Creation Steps

Create a shortcut with a custom icon to launch the application.

1.  **Create a new shortcut:** Right-click on your desktop or any desired location, then select "New" → **"Shortcut"**.
2.  **Enter the executable path:** In the "Type the location of the item" field, enter the full path to the executable and click "Next".
3.  **Name the shortcut:** Give the shortcut a name (e.g., `MoleditPy`) and click "Finish".
4.  **Change the icon:**
      * Right-click the newly created shortcut and select **"Properties"**.
      * In the "Shortcut" tab, click the **"Change Icon..."** button.
      * Browse to the icon file path mentioned above and select `icon.ico`.
      * Click "OK" → "Apply" → "OK" to finish.

### 4\. How to Pin the Shortcut to the Start Menu

To make the shortcut appear in the Windows Start Menu program list, place it in the following folder:

1.  **Locate the Start Menu folder:** Open File Explorer and enter one of the following paths in the address bar:
      * **For the current user only (Recommended):** `%APPDATA%\Microsoft\Windows\Start Menu\Programs`
      * **For all users (Requires administrator privileges):** `%ALLUSERSPROFILE%\Microsoft\Windows\Start Menu\Programs`
2.  **Copy or move the shortcut:** **Copy or move** the shortcut you created (e.g., the `MoleditPy` shortcut from your desktop) into this folder. (Holding `Ctrl` while dragging will copy it.)
3.  **Confirm:** Click the Start button and check that the shortcut name now appears in the program list.

-----

### Windows: ショートカットとファイル関連付けガイド

Windows環境でプロジェクトファイル（`.pmeraw`）をダブルクリックで開いたり、アイコン付きのショートカットからアプリケーションを起動したりするための詳細な設定手順です。

#### 1\. 実行ファイルとアイコンのパス

| 項目 | 詳細 |
| :--- | :--- |
| **実行ファイル (Target Application)** | `C:\Users\%USERNAME%\AppData\Local\Programs\Python\PythonXX\Scripts\moleditpy.exe` |
| **関連付けたい拡張子** | `.pmeraw` |
| **アイコンファイルのパス (Icon)** | `C:\Users\%USERNAME%\AppData\Local\Programs\Python\PythonXX\Lib\site-packages\moleditpy\assets\icon.ico` |

> **注意:** `%USERNAME%` はお使いのユーザー名、`PythonXX` はインストールされているPythonのバージョンディレクトリ名（例: `Python311`）に置き換えてください。

#### 2\. ファイルの関連付け手順 (.pmeraw をダブルクリックで開く)

`.pmeraw` ファイルを `moleditpy.exe` で開けるようにする手順です。

1.  **任意の `.pmeraw` ファイルを見つける:** エクスプローラーで、関連付けを行いたいファイル（例: `sample.pmeraw`）を**右クリック**します。
2.  **関連付けの変更:** コンテキストメニューから\*\*「プロパティ」**を選択し、「全般」タブの「ファイルの種類」の横にある**「変更...」\*\*ボタンをクリックします。
3.  **プログラムの選択:** 「このファイルを開く方法を選んでください」というウィンドウが表示されたら、「その他のアプリ」→\*\*「PCでアプリを探す」\*\*を選択します。
4.  **実行ファイルの指定:** 上記のパス（プレースホルダーを置き換えたもの）にある **`moleditpy.exe`** を指定し、「開く」をクリックして完了です。

> **補足:** カスタムアイコンをシステム全体のファイル関連付けに適用するには、通常、レジストリを編集する必要があります。レジストリ操作は上級者向けであり、慎重に行ってください。

#### 3\. デスクトップショートカットの作成手順

アプリケーションを起動するためのアイコン付きショートカットを作成します。

1.  **新規ショートカットの作成:** デスクトップなどの任意の場所を右クリックし、「新規作成」→\*\*「ショートカット」\*\*を選択します。
2.  **実行ファイルのパスを入力:** 「項目の場所を入力してください」に、上記の実行ファイルのパスを入力し、「次へ」をクリックします。
3.  **名前を指定:** ショートカットに任意の名前（例: `MoleditPy`）を付けて、「完了」をクリックします。
4.  **アイコンの変更:**
      * 作成されたショートカットを右クリックし、\*\*「プロパティ」\*\*を選択します。
      * 「ショートカット」タブの\*\*「アイコンの変更...」\*\*ボタンをクリックします。
      * 上記のアイコンファイルのパスを参照して `icon.ico` を選択します。
      * 「OK」→「適用」→「OK」をクリックして完了です。

#### 4\. ショートカットをスタートメニューに配置する方法

作成したショートカットをWindowsのスタートメニューのプログラム一覧に表示するには、以下の手順でフォルダに配置します。

1.  **スタートメニューフォルダのパスを確認:** エクスプローラーのアドレスバーに以下のいずれかのパスを入力し、Enterキーを押します。
      * **現在のユーザー専用 (推奨):** `%APPDATA%\Microsoft\Windows\Start Menu\Programs`
      * **すべてのユーザー共通 (管理者権限が必要):** `%ALLUSERSPROFILE%\Microsoft\Windows\Start Menu\Programs`
2.  **ショートカットのコピーまたは移動:** 上記で開いたスタートメニューフォルダに、作成したショートカット（例: デスクトップ上の `MoleditPy` ショートカット）を**コピーまたは移動**します。（Ctrlキーを押しながらドラッグするとコピーされます。）
3.  **確認:** スタートボタンをクリックし、プログラム一覧からショートカットに付けた名前が追加されていることを確認してください。

-----

## macOS: How to create a .app from a Python CLI app via Automator

This guide explains how to turn a Python app installed with `pip install moleditpy` into a native **.app** application on macOS, complete with a custom icon.

### 1\. Locate the Binary

The executable binary will be generated at a path similar to this:

```
/Users/<username>/Library/Python/<python_version>/bin/moleditpy
```

### 2\. Create a New Application in Automator

1.  Open **Automator**.
2.  Select **"New Document"** → **"Application"**.
3.  From the left pane, find "Utilities" → and add **"Run Shell Script"**.

### 3\. Configure the Script

Enter the following into the script area. The `"$@"` part is to pass any files dropped onto the app as arguments.

```bash
#!/bin/bash
/Users/<username>/Library/Python/<python_version>/bin/moleditpy "$@"
```

### 4\. Save as an Application

1.  From the menu, select **"File"** → **"Save"**.
2.  Name it **`MoleditPy.app`**.
3.  You can save it to your Desktop temporarily.

### 5\. Move to the Applications Folder

After saving, copy the created app to the system's standard Applications folder.

**Copying via Finder:**

1.  Select `MoleditPy.app` in Finder and copy it with `⌘ + C`.
2.  Open the **"Applications"** folder from the Finder's "Go" menu.
3.  Paste it with `⌘ + V`.

**Copying via Terminal:**

```bash
sudo cp -R ~/Desktop/MoleditPy.app /Applications/
```

### 6\. Set the Icon

Location of the icon image:

```
/Users/<username>/Library/Python/<python_version>/lib/python<python_version>/site-packages/moleditpy/assets/icon.png
```

**Icon Setting Steps (Copy-Paste Method):**

1.  Open the `icon.png` from the path above in Finder.
2.  With the image open in "Preview," press `⌘ + A` to **select all**, then `⌘ + C` to **copy**.
3.  In Finder, select `MoleditPy.app` in the `/Applications` folder and press `⌘ + I` (Get Info).
4.  Click the small icon in the top-left corner (a border will appear) → press `⌘ + V` to **paste**.
    The app's icon will now be updated.

### 7\. Check Execution Permissions (if necessary)

If the app fails to launch, run the following command:

```bash
chmod +x /Users/<username>/Library/Python/<python_version>/bin/moleditpy
```

### 8\. Done\!

Now you can run the Python app just by double-clicking **`MoleditPy.app`** installed in `/Applications`.

-----

### macOS: Python CLIアプリをAutomator経由で.app化する方法

この手順では、`pip install moleditpy` でインストールしたPythonアプリをmacOS上で\*\*.appアプリケーション\*\*として使えるようにし、独自アイコンを設定します。

#### 1\. バイナリの場所を確認

以下のようなパスに実行ファイル（バイナリ）が生成されます：

```
/Users/<username>/Library/Python/<python_version>/bin/moleditpy
```

#### 2\. Automatorで新規アプリケーションを作成

1.  **Automator** を開きます。
2.  **「新規書類」** → **「アプリケーション」** を選択します。
3.  左の一覧から「ユーティリティ」→\*\*「シェルスクリプトを実行」\*\* を追加します。

#### 3\. スクリプトを設定

スクリプト欄に以下を入力します。`"$@"` は、ドラッグ＆ドロップされたファイルを引数として渡すためのものです。

```bash
#!/bin/bash
/Users/<username>/Library/Python/<python_version>/bin/moleditpy "$@"
```

#### 4\. アプリとして保存

1.  メニューから\*\*「ファイル」\*\* → **「保存」** を選択します。
2.  名前を **`MoleditPy.app`** にします。
3.  保存場所は一時的にデスクトップなど任意で構いません。

#### 5\. アプリケーションフォルダーへ移動

保存後、作成したアプリをシステム標準のアプリケーションフォルダーへコピーします。

**Finderからコピーする場合**

1.  Finderで `MoleditPy.app` を選択し、`⌘ + C` でコピーします。
2.  Finderメニューの「移動」→\*\*「アプリケーション」\*\* を開きます。
3.  `⌘ + V` で貼り付けます。

**ターミナルからコピーする場合**

```bash
sudo cp -R ~/Desktop/MoleditPy.app /Applications/
```

#### 6\. アイコンを設定

アイコン画像の場所：

```
/Users/<username>/Library/Python/<python_version>/lib/python<python_version>/site-packages/moleditpy/assets/icon.png
```

**アイコン設定手順（コピペ方式）**

1.  Finderで上記の `icon.png` を開きます。
2.  「プレビュー」で画像を開いた状態で `⌘ + A` で**全選択**、`⌘ + C` で**コピー**します。
3.  Finderで `/Applications` フォルダー内の `MoleditPy.app` を選択 → `⌘ + I`（情報を見る）を押します。
4.  左上の小さいアイコンをクリック（枠が出る） → `⌘ + V` で**貼り付け**ます。
    これでアプリのアイコンが変更されます。

#### 7\. 実行権限を確認（必要な場合）

起動できない場合は、以下を実行してください。

```bash
chmod +x /Users/<username>/Library/Python/<python_version>/bin/moleditpy
```

#### 8\. 完了

これで `/Applications` にインストールされた **`MoleditPy.app`** をダブルクリックするだけでPythonアプリが実行できます。

-----

## Technical Details

  * **GUI and 2D Drawing (PyQt6):**
      * Interactively manipulates custom `AtomItem` and `BondItem` objects on a `QGraphicsScene`.
      * The Undo/Redo feature is implemented by serializing the entire application state with `pickle` and storing it on a stack.
  * **Chemical Calculations (RDKit / Open Babel):**
      * Generates RDKit molecule objects from 2D data to perform 3D coordinate generation (`AllChem.EmbedMolecule`) and calculate molecular properties.
      * When 3D coordinate generation with RDKit fails, it falls back to **Open Babel** to attempt the calculation.
      * Calculations are performed on a separate thread (`QThread`) to maintain GUI responsiveness.
  * **3D Visualization (PyVista / pyvistaqt):**
      * Generates and renders PyVista meshes (spheres and cylinders) from RDKit conformer coordinates.
      * Implements a custom `vtkInteractorStyle` to enable direct drag-and-drop editing of atoms within the 3D view.

### 技術的な仕組み

  * **GUI と 2D 描画 (PyQt6):**
      * `QGraphicsScene` 上にカスタムの `AtomItem`（原子）と `BondItem`（結合）を配置し、対話的に操作します。
      * Undo/Redo機能は、アプリケーションの状態を丸ごと `pickle` でシリアライズしてスタックに保存することで実現しています。
  * **化学計算 (RDKit / Open Babel):**
      * 2D データから RDKit 分子オブジェクトを生成し、3D 座標生成（`AllChem.EmbedMolecule`）や分子特性計算を実行します。
      * RDKitでの3D座標生成が失敗した際は、**Open Babel**にフォールバックして計算を試みます。
      * 計算は別スレッド（`QThread`）で行い、GUI の応答性を維持しています。
  * **3D 可視化 (PyVista / pyvistaqt):**
      * RDKit のコンフォーマ座標から PyVista のメッシュ（球や円柱）を生成して描画します。
      * カスタムの`vtkInteractorStyle`を実装し、3Dビュー内での原子の直接的なドラッグ＆ドロップ編集を可能にしています。

-----

## License

This project is licensed under the **Apache-2.0 License**. See the `LICENSE` file for details.


### ライセンス

このプロジェクトは **Apache-2.0 License** のもとで公開されています。詳細は `LICENSE` ファイルを参照してください。
