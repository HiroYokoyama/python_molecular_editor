# Linux環境でmoleditpyアプリを安定動作させるための開発検討（まだ動きません）

PyQt6, PyVista, RDKit, Open Babelなど、C++ライブラリに依存する複雑なPythonスタックをLinux環境で安定して動作させるための、推奨される環境構築手順です。

`pip`によるインストールは、ライブラリが独自に同梱する共有ライブラリ間の競合（`Segmentation fault`の原因）を引き起こす可能性があるため、以下の`conda`を用いた方法を強く推奨します。

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
