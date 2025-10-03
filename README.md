# moleditpy -- Python Molecular Editor

Pythonで構築された、シンプルで直感的な分子構造エディターです。2Dでの分子描画と、3D構造可視化をサポートします。

作者: HiroYokoyama
ライセンス: Apache-2.0
リポジトリ: [https://github.com/HiroYokoyama/python\_molecular\_editor](https://github.com/HiroYokoyama/python_molecular_editor)

![](img/screenshot.png)

-----

## 概要

このアプリケーションは、化学者や学生が分子構造を容易に描き、その3次元的な形状を視覚的に確認するためのツールです。PyQt6によるモダンなGUI、RDKitによる強力な化学計算、PyVistaによる高性能な3Dレンダリングを組み合わせています。

-----

## 主な機能

  * **直感的な2D描画**

      * マウスのクリック＆ドラッグで原子や結合を簡単に追加・編集
      * 原子（C, N, O, Hなど）や結合（単結合、二重結合、三重結合）のツールバーからの選択
      * 周期表ダイアログから任意の元素を選択可能
      * Undo/Redo、選択、全選択、削除

  * **テンプレートプレビューと配置**

      * ベンゼン環や3〜9員環のテンプレートをプレビューして配置可能
      * 既存の原子や結合にスナップして配置できる

  * **キーボードショートカットと操作性の改善**

      * `Space`: 選択モード切替 / 選択モードで全選択
      * `1`/`2`/`3`: カーソル下または選択中の結合の結合次数を変更
      * `Delete` / `Backspace`: 選択項目の削除
      * 原子上でキー入力（`C`, `N`, `O`, `S`, `F`, `B`, `I`, `H`, `Shift+C`=Cl, `Shift+B`=Br）で元素を即時切替

  * **2D構造の最適化**

      * RDKit の `Compute2DCoords` を使った自動レイアウト（Optimize 2D）

  * **高品質な3D可視化**

      * RDKit で 3D 座標を生成し MMFF94 ベースで最適化
      * PyVista / pyvistaqt によるインタラクティブな3D表示（原子は球、結合は円柱で表現）
      * 3Dビューをクリックすると2Dにフォーカスを戻す動作など、操作性の改善

  * **ファイル入出力**

      * 2D を MOL 形式で保存
      * 3D を MOL / XYZ 形式で保存
      * MOL/SDF の読み込み
      * プロジェクトファイル（`.pmeraw`）で編集状態の保存/読み込み

-----

## 実行とインストール

#### 必要ライブラリ

`PyQt6`, `RDKit`, `NumPy`, `PyVista`, `pyvistaqt`

#### インストール例

**pip を使う場合:**

```bash
pip install moleditpy
```

> **Note**
> RDKit は `conda` を使ってインストールすることが推奨されます。

#### アプリの起動

```bash
moleditpy
```

-----

## 技術的な仕組み

  * **GUI と 2D 描画 (PyQt6)**

      * アプリの骨格と 2D 描画キャンバスは PyQt6 の Graphics View Framework を利用しています。
      * `QGraphicsScene` 上にカスタムの `AtomItem`（原子）と `BondItem`（結合）を配置し、マウス/キー入力で対話的に操作します。
      * テンプレートのプレビューは専用の `TemplatePreviewItem` で描画されます。

  * **化学計算 (RDKit)**

      * 2D 描画データから MOL ブロックを生成し、RDKit に渡して 3D 座標生成（`AllChem.EmbedMolecule`）と最適化（`AllChem.MMFFOptimizeMolecule`）を実行します。
      * 2D の自動レイアウト（`Compute2DCoords`）にも RDKit を使用します。
      * 計算は別スレッド（`QThread`）で行い、GUI の応答性を維持しています。

  * **3D 可視化 (PyVista / pyvistaqt)**

      * RDKit のコンフォーマ座標から PyVista のメッシュ（球や円柱）を生成して描画します。
      * ボンドの種類（単結合/二重結合/三重結合）に応じた複数円柱のレンダリングを実装しています。

-----

## ライセンス

このプロジェクトは Apache-2.0 License のもとで公開されています。詳細は `LICENSE` ファイルを参照してください。
