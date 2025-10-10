# Linux向けインストール

```bash
pip install moleditpy-linux
```
ライブラリの競合により、Segmentation faultが起こるので、OpenBabelを無効にしたバージョン。

### 具体的な方法

```bash
python3 -m venv moleditpy
source moleditpy/bin/activate
pip install moleditpy-linux
sudo apt install libxcb-cursor0
moleditpy
```

-----

# 検討メモ

### Issue: Linux環境における起動時クラッシュ問題の調査と対応

#### 1. 問題の概要

特定のLinux環境において、アプリケーション `moleditpy` を起動しようとすると、メインウィンドウが表示される前にクラッシュする。

- **発生環境**: 特定のLinux環境
- **正常動作環境**: Windows, macOSでは問題なく起動することが確認されている。
- **エラーメッセージ**: `X Error of failed request: BadWindow (invalid Window parameter)`

#### 2. 調査の経緯

`print`文を用いたトレースにより、`MainWindow`のUI初期化処理中、**`pyvista.QtInteractor` ウィジェットを生成する行でクラッシュが発生すること**を特定した。

#### 3. 原因の特定

直接的な原因は、**PyVistaのQtウィジェットである `QtInteractor` の初期化失敗**にある。

Windows/Macで動作することから、アプリケーションのコードロジックではなく、Linux環境固有の以下の要因が複合的に絡んでいる問題と推測される。

- **グラフィックドライバ**と **OpenGL** の互換性
- **PyQt6** と **VTK** (PyVistaのバックエンド) のライブラリバージョンの競合
- Linuxの **X11ウィンドウシステム** との相性

#### 4. 対応と現状

現状、根本的な解決には至っておらず、複数の問題点を回避するために以下の検討を行った。

- **Open Babel機能の無効化**:
  **`openbabel` (pybel) をインポートすると、UI初期化以前のimport時点でアプリケーションがクラッシュする**。このため、`import`文を含め、`openbabel`関連の機能はすべてコメントアウトして検討。

- **仮想環境 (`venv`) での解決試行**:
  依存関係の競合を解決するため、Python仮想環境 (`venv`) を構築し、クリーンな状態でライブラリを再インストールしたが、**問題は解決しなかった**。これにより、問題がPythonライブラリのバージョンだけでなく、より低レイヤーのシステムライブラリやドライバに起因する可能性が非常に高いことが示唆される。

#### 5. 今後のタスク

環境依存の問題である可能性が非常に高いため、以下のシステムレベルでの調査が必要となる。


-----


## `moleditpy` のインストール時に発生したQtプラットフォームプラグインエラーのメモ

## 1\. 発生した問題

Condaで新しい環境を作成し、`moleditpy`をインストールして実行しようとしたところ、GUIが起動せずエラーが発生した。

### 1.1. 初期インストール手順

以下のコマンドで環境を構築した。

```bash
# 1. 新しい環境を作成
conda create -n moledit_env python=3.11 -c conda-forge

# 2. 環境を有効化
conda activate moledit_env

# 3. 必要なライブラリをCondaでインストール
conda install -c conda-forge pyqt rdkit numpy pyvista pyvistaqt openbabel

# 4. moleditpyをpipでインストール
pip install moleditpy

# 5. 実行
moleditpy
```

### 1.2. 最初に発生したエラー

```
qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
```

### 1.3. 発生したエラー詳細（QT_DEBUG_PLUGINS=1）


```
/.../libQt6XcbQpa.so.6: undefined symbol: _ZN20QSpiAccessibleBridgeC1Ev, version Qt_6_PRIVATE_API
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized.
```

## 2\. 原因の分析

1.  **初期エラー (`libxcb-cursor0 is needed`)**: QtがLinux上でGUIを描画するために必要な、OSレベルの共有ライブラリが不足していることが原因。
2.  **二次エラー (`undefined symbol`)**: より深刻な問題を示唆。これはOSライブラリの不足ではなく、**Conda環境内にインストールされたQt関連ライブラリ同士のバージョンやビルドに不整合が生じている**ことが原因。`conda`でインストールした`pyqt`と、`pip`経由で構築された環境との間でABI（Application Binary Interface）の互換性が崩れた可能性が高い。


PyQt6, PyVista, RDKit, Open Babelなど、C++ライブラリに依存する複雑なPythonスタックをLinux環境で安定して動作させるための、推奨される環境構築手順です。

`pip`によるインストールは、ライブラリが独自に同梱する共有ライブラリ間の競合（`Segmentation fault`の原因）を引き起こす可能性があるため、以下の`conda`を用いた方法を強く推奨します。


-----

## メモ2

### 前提条件

  - `miniconda` または `anaconda` がインストールされていること。

### 手順

1.  **新しいConda仮想環境を作成する**
    ターミナルを開き、以下のコマンドでアプリケーション専用の独立した環境を作成します。Pythonバージョンは3.11など、安定したLTS版を推奨します。

    ```bash
    # "moledit-env" という名前でPython 3.11の環境を新規作成
    conda create -n moledit-env python=3.11 -y
    ```

2.  **作成した環境を有効化する**

    ```bash
    conda activate moledit-env
    ```

    プロンプトの先頭が `(moledit-env)` のように変わります。

3.  **全ての依存ライブラリを`conda-forge`から一度にインストールする**
    必要なライブラリを一つのコマンドにまとめてインストールすることで、`conda`が互換性のあるバージョンを解決してくれます。

    ```bash
    conda install -c conda-forge pyqt pyvista vtk rdkit openbabel -y
    ```

4.  **アプリケーションを実行する**
    全ての準備が整いました。このクリーンな環境でアプリケーションを実行します。

    ```bash
    sudo apt install libxcb-cursor0
    pip install moleditpy
    moleditpy
    ```

上記の手順に従うことで、本ドキュメント下部に詳述するような共有ライブラリの競合問題を未然に防ぎ、安定した動作が期待できます。

-----

## 詳細: Linux環境におけるSegmentation Faultの調査経緯と原因

### 要約

PyQt6, PyVista, RDKit, Open Babel を使用したPython GUIアプリケーションが、Linux環境で起動時に `Segmentation fault` でクラッシュした。調査の結果、原因は **`pip`でインストールされた`openbabel-wheel`が、PyQt6と互換性のない共有ライブラリ (`libxcb.so`) を独自に同梱しており**、ライブラリのバージョン競合を引き起こしていたことであった。

**解決策**として、`pip`の使用をやめ、**`conda` (`conda-forge`チャンネル) を使って全ての依存ライブラリをクリーンな仮想環境にインストールする**ことで、ライブラリの依存関係が正しく管理され、問題が解決した。

### 現象

開発した分子エディタ (`main.py`) をLinux環境で実行したところ、ウィンドウが表示されることなく、以下のエラーで即座にプロセスが終了する。

```bash
$ python3 main.py
Segmentation fault (core dumped)
```

### 環境

  - **OS**: Linux
  - **主なライブラリ**:
      - Python 3.10+
      - PyQt6 (GUIフレームワーク)
      - PyVista, VTK (3D描画)
      - RDKit, Open Babel (化学計算)
  - **パッケージ管理**: miniconda3, pip

### 調査の経緯

`Segmentation fault` は低レベルなメモリエラーであり、Pythonコードのロジックよりも、C++で書かれたライブラリ間の競合や環境要因が原因であることが多い。以下の手順で原因を絞り込んだ。

1.  **環境要因の切り分け**: `Wayland` vs `X11` の問題や、`libxcb-cursor0` のような必須ライブラリの不足を疑ったが、いずれも直接の原因ではなかった。

2.  **GDBによるバックトレースの取得**: GDB (GNUデバッガ) を使用してクラッシュした瞬間の詳細なスタックトレースを取得した。

    ```c
    #0  0x00007fffbec0de05 in _xcb_in_expect_reply () from /home/user/miniconda3/lib/python3.13/site-packages/openbabel/lib/openbabel/3.1.0/../../../../openbabel_wheel.libs/libxcb-65da195c.so.1.1.0
    ...
    #10 0x00007fffc95e817a in QGuiApplicationPrivate::createPlatformIntegration() () from /home/user/miniconda3/lib/python3.13/site-packages/PyQt6/Qt6/lib/libQt6Gui.so.6
    ...
    #14 0x00007fffca198a6d in QApplicationPrivate::init() () from /home/user/miniconda3/lib/python3.13/site-packages/PyQt6/Qt6/lib/libQt6Widgets.so.6
    ```

### 根本原因

GDBのバックトレース（特にフレーム `#0`）から、決定的な証拠が判明した。`QApplication` の初期化プロセス中に呼び出された `libxcb.so` が、システムやcondaの標準ライブラリではなく、**`openbabel` の `pip` wheelパッケージ内部に同梱されていたもの**だった。

この `openbabel` 由来の `libxcb.so` が、PyQt6が期待するバージョンと異なっていたため、ライブラリ間で深刻な競合が発生し、`Segmentation fault` を引き起こしていた。

### 結論

調査の結果、本問題の唯一の確実な解決策は、ドキュメント冒頭に記載した**condaを用いた環境の再構築**であることが判明した。

科学技術計算系のPythonライブラリ（特にGUIや3D描画が絡むもの）をLinux環境で使用する際、`pip`でインストールされたwheelパッケージがシステムライブラリと競合することがある。このような複雑な依存関係を持つスタックでは、`conda` (`conda-forge`チャンネル) を用いて環境を構築することが、安定した動作を実現するための最も確実な方法である。
