# MoleditPy 5 への移行

MoleditPy 5.0.0a1 は Windows・macOS・Linux で共通の `moleditpy` パッケージを使用します。Python 3.9–3.14、PySide6-Essentials、RDKit、PyVista が必要です。コアは Open Babel をインストール・読み込み・バンドルしません。個別の Linux パッケージとビルド・公開処理は廃止しました。

```sh
python -m pip uninstall moleditpy-linux openbabel-wheel
python -m pip install --upgrade --pre 'moleditpy>=5.0.0a1,<6'
```

PyQt6 からの移行は新しい仮想環境を推奨します。同一プロセスで Qt バインディングを混在させないでください。Open Babel は他のプラグインやアプリが必要としていなければ削除できます。

Fallback は RDKit の後に Direct を試します。Direct は描画した座標を維持して水素を追加します。保存済み Open Babel 変換設定は Fallback へ、旧組み込み Open Babel 最適化設定は RDKit MMFF94s へ移行します。追加の力場は、既存の登録 API を使用する任意の最適化プラグインで提供できます。

プラグインは PySide6 の `Signal`・`Slot`・`Property` と Shiboken6 の有効性確認へ移行してください。ワーカースレッドから UI を更新する場合は GUI スレッドに属する QObject と型付きスロットを使用します。

コアは LGPL-3.0-or-later です。独立したプラグインは任意のライセンスを選択でき、既存プラグインのライセンスは変更しません。詳細は [LICENSING.md](../LICENSING.md) と [DISTRIBUTION.md](../DISTRIBUTION.md) を参照してください。このブランチはリリースや PyPI 公開を行いません。
