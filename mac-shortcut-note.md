# macOSでPython CLIアプリをAutomator経由で.app化する方法（moleditpy）

この手順では、`pip install moleditpy` でインストールしたPythonアプリをmacOS上で`.app`アプリケーションとして使えるようにし、独自アイコンを設定します。

---

## 1. バイナリの場所を確認

以下のようなパスに実行ファイル（バイナリ）が生成されます：

```

/Users/<username>/Library/Python/<python_version>/bin/moleditpy

````

---

## 2. Automatorで新規アプリケーションを作成

1. **Automator** を開く  
2. 「新規書類」→「アプリケーション」を選択  
3. 左の一覧から「ユーティリティ」→「シェルスクリプトを実行」を追加  

---

## 3. スクリプトを設定

スクリプト欄に以下を入力します：

```bash
#!/bin/bash
/Users/<username>/Library/Python/<python_version>/bin/moleditpy "$@"
````

> `"$@"` は、ドラッグ＆ドロップされたファイルなどを引数として渡すためのものです。

---

## 4. アプリとして保存

1. メニューから「ファイル」→「保存」
2. 名前を `moleditpy.app` にする
3. 保存場所は一時的にデスクトップなど任意でOKです

---

## 5. アプリケーションフォルダーへ移動

保存後、作成したアプリをシステム標準のアプリケーションフォルダーへコピーします。

### Finderからコピーする場合

1. Finderで `moleditpy.app` を選択
2. `⌘ + C` でコピー
3. Finderメニューの「移動」→「アプリケーション」を開く
4. `⌘ + V` で貼り付け

### ターミナルからコピーする場合

```bash
sudo cp -R ~/Desktop/moleditpy.app /Applications/
```

> `sudo` 実行時にパスワード入力を求められる場合があります。

---

## 6. アイコンを設定

アイコン画像の場所：

```
/Users/<username>/Library/Python/<python_version>/lib/python<python_version>/site-packages/moleditpy/assets/icon.png
```

### アイコン設定手順（コピペ方式）

1. Finderで上記の`icon.png`を開く
2. 「プレビュー」で画像を開いた状態で

   * `⌘ + A` で全選択
   * `⌘ + C` でコピー
3. Finderで `/Applications` フォルダー内の `moleditpy.app` を選択 → `⌘ + I`（情報を見る）
4. 左上の小さいアイコンをクリック（枠が出る） → `⌘ + V` で貼り付け

これでアプリのアイコンが変更されます。


### PNGをICNS形式に変換して設定する場合：
```
# PNG → ICNS 変換
iconutil -c icns /path/to/icon.iconset

# アプリに適用
cp /path/to/icon.icns /Applications/MyTool.app/Contents/Resources/app.icns
備考: iconutil は macOS 標準コマンドです。icon.iconset は以下のように作成できます：
mkdir icon.iconset
cp icon.png icon.iconset/icon_512x512.png
iconutil -c icns icon.iconset
```
---

## 7. 実行権限を確認（必要な場合）

起動できない場合は、以下を実行：

```bash
chmod +x /Users/<username>/Library/Python/<python_version>/bin/moleditpy
```

---

## 8. 完了

これで `/Applications` にインストールされた `moleditpy.app` をダブルクリックするだけでPythonアプリが実行できます。


---
