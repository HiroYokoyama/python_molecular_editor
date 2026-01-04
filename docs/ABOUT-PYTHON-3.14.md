# Python Version Support / Python バージョンのサポートについて

> [!NOTE]
> Python 3.14 Support Status
> Python 3.14 のサポート状況について

## English

### Compatibility Notice
Moleditpy currently fully supports Python 3.9 through 3.13.
Support for Python 3.14 is currently pending.

### Technical Details
This application relies on openbabel-wheel for core functionality.
Currently, pre-compiled binary wheels compatible with Python 3.14 are not yet available from the upstream source.

Installing on Python 3.14 will trigger a compilation from source, which may fail due to complex system requirements. To ensure a stable experience, installation is restricted to Python 3.13 and earlier until compatible wheels become available.

### Recommendation
Please use Python 3.13 (or earlier).
We will update our dependency requirements as soon as the upstream ecosystem supports Python 3.14.

---

## 日本語

### 互換性のお知らせ
Moleditpy は現在、Python 3.9 から 3.13 までの動作を完全にサポートしています。
Python 3.14 への対応は、現在待機中です。

### 技術的な背景
本アプリケーションは openbabel-wheel に依存していますが、現時点では Python 3.14 と互換性のあるバイナリ（Wheel）がまだ利用可能になっていません。

Python 3.14 環境では、ソースコードからのビルドが必要となり、インストールが不安定になる可能性があります。そのため、互換性のあるバイナリが公開されるまでの間、安全のために Python 3.13 以下への制限を設けています。

### 推奨環境
Python 3.13（またはそれ以前）をご利用ください。
上流のライブラリ環境が Python 3.14 に対応次第、速やかに制限を解除する予定です。
