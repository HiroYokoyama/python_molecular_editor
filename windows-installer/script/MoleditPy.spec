# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path


spec_dir = Path(SPECPATH).resolve()
repo_root = spec_dir.parent.parent
linux_package_dir = repo_root / 'moleditpy-linux' / 'src' / 'moleditpy_linux'

a = Analysis(
    [str(linux_package_dir / '__main__.py')],
    pathex=[str(linux_package_dir)],
    binaries=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
    datas=[
        (str(linux_package_dir / 'assets'), 'assets'),
    ],
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='MoleditPy',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=[str(linux_package_dir / 'assets' / 'icon.ico')],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='MoleditPy',
)
