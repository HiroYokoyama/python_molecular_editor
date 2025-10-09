### English
# ⚠️ Notice Regarding Python 3.14 Support
This project currently does not support Python 3.14.
Users attempting to install or run this project with Python 3.14 will likely encounter installation failures.
The Reason
The core issue lies with the broader scientific Python ecosystem. Major dependencies, such as NumPy, and others with C-extensions, have not yet released stable, pre-compiled binary packages (wheels) for Python 3.14.
When pip cannot find a compatible wheel, it falls back to building the package from source. This process requires a C/C++ compiler and a full development environment (e.g., Visual Studio Build Tools on Windows), which most users do not have configured. This results in compilation errors, such as Unknown compiler(s): [['cl'], ['cc'], ...].
Recommendations
For Users:
We strongly recommend using a stable, fully supported version of Python. As of now, Python 3.12 is the recommended version for this project to ensure a smooth installation experience. Please use a virtual environment to manage your Python versions.
For Contributors:
To prevent users from running into this issue, it is best practice to explicitly define the supported Python versions in your pyproject.toml:

I will monitor the status of the ecosystem and begin official support for Python 3.14 once the core dependencies are fully available and stable.


### 日本語 (Japanese)
# ⚠️ Python 3.14 のサポートに関する注意喚起
現在、このプロジェクトは Python 3.14 をサポートしていません。
Python 3.14 を使用してこのプロジェクトのインストールや実行を試みると、インストールの失敗に繋がる可能性が非常に高いです。
理由
根本的な問題は、Pythonの科学技術計算エコシステム全体の対応遅延にあります。NumPy, SciPy, Pandasといった、C言語の拡張モジュールを含む主要な依存ライブラリが、Python 3.14に対応した安定版のコンパイル済みバイナリパッケージ（Wheel）をまだリリースしていません。
pipは対応するWheelを見つけられない場合、ソースコードからパッケージをビルドしようとします。このプロセスにはC/C++コンパイラと開発環境（例：WindowsにおけるVisual Studio Build Tools）が必要であり、ほとんどのユーザー環境には設定されていません。結果として、Unknown compiler(s): [['cl'], ['cc'], ...] のようなコンパイルエラーが発生します。
さらに、複雑な依存関係を解決しようとする過程で、pipが制約を満たすために古いバージョンのパッケージをインストールしようとし、それも同じ理由でビルドに失敗するという事態も起こり得ます。
推奨事項
ユーザーの方へ：
安定しており、完全に対応が済んでいるPythonバージョンを使用することを強く推奨します。現時点では、スムーズなインストールが可能な Python 3.13 がこのプロジェクトの推奨バージョンです。仮想環境を使ってPythonのバージョンを管理してください。
コントリビューター（開発協力者）の方へ：
ユーザーがこの問題に直面するのを防ぐため、pyproject.toml に対応するPythonバージョンを明記することがベストプラクティスです。

エコシステムの状況を注視し、主要な依存関係が安定し利用可能になり次第、Python 3.14の公式サポートを開始する予定です。
