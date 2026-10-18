# Migrating to MoleditPy 5

MoleditPy 5 uses one `moleditpy` package on Windows, macOS, and Linux. It requires
Python 3.10–3.14 and official `openbabel >= 3.2.1, < 3.3` wheels. The separate
`moleditpy-linux` distribution is retired.

In the Python environment containing the old application, run:

```sh
python -m pip uninstall moleditpy-linux openbabel-wheel
python -m pip install --upgrade 'moleditpy>=5,<6'
```

Removing the old distributions first avoids overlapping `openbabel` files and
the shared `moleditpy` command. Start the application with `moleditpy`, or use
`python -m moleditpy_installer` to recreate shortcuts and file associations.
The installer helper remains a shortcut creator; the frozen application builds
come from the main `moleditpy` package.

The official Open Babel wheels cover Linux x86-64 and ARM64 with glibc 2.17+.
The full application also needs compatible Qt, RDKit, and VTK libraries; the
Open Babel glibc minimum alone does not establish support for a Linux system.
Alpine/musl, 32-bit Linux, free-threaded Python, and Python 3.15 are outside the
tested support matrix.

Open Babel loads native symbols globally on Linux. Load RDKit before Open Babel
in plugins and scripts that use both. MoleditPy follows that order and checks
conversion, bundled force-field parameters, and cold GUI startup in CI.

Sources: [official Open Babel wheels](https://pypi.org/project/openbabel/),
[wheel deprecation notice](https://github.com/njzjz/openbabel-wheel),
[native import-order troubleshooting](https://gentle1999.github.io/MolOP/0.2.13/en/guides/troubleshooting/).

## Release and PyPI checklist

1. Publish `MoleditPy 5.0.0` from the unified release workflow after CI passes.
2. Publish the updated shortcut helper, which depends on `moleditpy >= 5, < 6`
   on every supported OS.
3. Archive the old `MoleditPy-linux` project from its PyPI settings. Keep all
   existing releases available for users pinned to older versions. Archival
   preserves downloads and prevents new uploads; it does not migrate existing
   installations automatically.

The existing `MoleditPy` PyPI project is reused. Its Trusted Publisher remains
the release workflow; there is no new package name or publisher to register.
See [PyPI project archival](https://blog.pypi.org/posts/2025-01-30-archival/).
