### Condaを使ったGUIソフト`moleditpy`の環境構築・実行手順書

この手順書は、GUIソフトウェアである`moleditpy`を、パッケージマネージャー`Conda`を用いて安定した実行環境にインストールし、起動するためのものです。

開発環境と必要なライブラリのバージョンを可能な限り再現し、GUIアプリケーションで特に問題となりやすいライブラリ間の互換性を重視した手順となっています。

#### 前提条件

  * Anaconda または Miniconda がPCにインストールされていること。

-----

### ステップ1: Conda環境の新規作成と有効化

まず、`moleditpy`専用のクリーンな仮想環境を作成します。

1.  **Condaプロンプトを開く**
    Windowsのスタートメニューから「Anaconda Prompt」または「Miniconda Prompt」を探して起動します。

2.  **新しい環境を作成する**
    互換性が高い**Python 3.12**を指定して、`moleditpy-env`という名前の環境を作成します。

    ```bash
    conda create -n moleditpy-env -c conda-forge python pyqt=6 rdkit openbabel pyvista matplotlib numpy pyvistaqt "vtk=9.4"
    conda install -c conda-forge pyqt=6
    ```

      * 実行中に `Proceed ([y]/n)?` と聞かれたら、`y` を入力してEnterキーを押します。

3.  **作成した環境を有効化する**
    これからの作業は、すべてこの新しい環境内で行います。

    ```bash
    conda activate moleditpy-env
    ```

    コマンドプロンプトの行の先頭が `(moleditpy-env)` に変わったことを確認してください。

-----


### ステップ2: `moleditpy`のインストール

必要な依存関係がすべて整いましたので、最後に`moleditpy`本体を`pip`でインストールします。

1.  **pipで`moleditpy`をインストール**
    `conda`で有効化した環境の中で、以下のコマンドを実行します。

    ```bash
    pip install moleditpy
    ```

-----

### ステップ3: `moleditpy`の起動

インストールが完了したら、`moleditpy`をGUIアプリケーションとして起動します。

1.  **起動コマンドの実行**
    コマンドプロンプト（`moleditpy-env`が有効化された状態）で以下を実行してみてください。

    ```bash
    moleditpy
    ```

    `moleditpy`のGUIウィンドウが立ち上がれば、環境構築はすべて成功です。


-----

### トラブルシューティング

  * **コマンドが見つからない(`command not found`)と表示された場合**: ステップ3で`moleditpy`のインストールが正しく完了しているか確認してください。
  * **起動時にエラーが出る場合**: 表示されたエラーメッセージを元に、どのライブラリで問題が起きているかを確認します。Conda環境を一度削除 (`conda env remove -n moleditpy-env`) して、手順を最初からやり直すと解決することも多いです。

今後は、`moleditpy`を使用する際に、必ず `conda activate moleditpy-env` コマンドでこの環境を有効化してから起動コマンドを実行してください。
