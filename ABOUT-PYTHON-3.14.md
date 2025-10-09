English
# ⚠️ Notice Regarding Python 3.14 Installation Issues
Users attempting to install this project with Python 3.14 may encounter installation failures.
## The Reason
The issue stems from a specific dependency conflict. The openbabel-wheel package, a dependency of this project, requires an older version of NumPy (e.g., numpy<2.0).
On Python 3.14, there are no pre-compiled binary packages (wheels) available for these older NumPy versions.
Consequently, pip attempts to build this older NumPy from source. This process requires a C/C++ compiler (like Visual Studio Build Tools on Windows) and will fail on most standard user environments, leading to a compilation error.
## Recommendation
To avoid this compilation issue, we strongly recommend using Python 3.13 for this project.
With Python 3.13, pip can find pre-compiled wheels for both openbabel-wheel and the specific older version of NumPy it requires, ensuring a smooth installation without the need for a local compiler environment.


日本語 (Japanese)
# ⚠️ Python 3.14 でのインストールに関する注意喚起
Python 3.14 環境でこのプロジェクトをインストールしようとすると、問題が発生する可能性があります。
## 理由
この問題は、特定の依存関係の競合が原因です。当プロジェクトが依存しているopenbabel-wheelパッケージが、古いバージョンのNumPy（例：numpy<2.0）を要求します。
しかし、Python 3.14環境では、この古いバージョンのNumPyに対応したコンパイル済みのバイナリパッケージ（Wheel）が提供されていません。
その結果、pipは古いNumPyをソースコードからビルドしようと試みます。このプロセスにはC/C++コンパイラ（WindowsにおけるVisual Studio Build Toolsなど）が必要であり、多くのユーザー環境ではこのビルドが失敗し、コンパイルエラーに至ります。
## 推奨事項
このコンパイル問題を回避するため、このプロジェクトでは Python 3.13 を使用することを強く推奨します。
Python 3.13環境であれば、pipがopenbabel-wheelと、それが必要とする古いNumPyの両方に対応したコンパイル済みWheelを見つけることができるため、ローカルでのコンパイル作業なしにスムーズなインストールが可能です。
