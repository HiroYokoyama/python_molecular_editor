# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path


spec_dir = Path(SPECPATH).resolve()
repo_root = spec_dir.parent.parent
linux_package_dir = repo_root / 'moleditpy-linux' / 'src' / 'moleditpy_linux'
src_dir = linux_package_dir.parent  # moleditpy-linux/src/ — makes `import moleditpy_linux` work

a = Analysis(
    [str(spec_dir / '_launcher.py')],
    pathex=[str(src_dir)],
    binaries=[],
    hiddenimports=[
        # rdkit.Chem.rdForceFieldHelpers is a C extension and imports these at
        # runtime, where the module graph cannot see it. Nothing in the reachable
        # Python graph pulls them in since 4.8.3 stopped importing AllChem, so the
        # frozen app died on its first import without them.
        'rdkit.ForceField',
        'rdkit.ForceField.rdForceField',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
    datas=[
        (str(linux_package_dir / 'assets'), 'assets'),
        (str(linux_package_dir / 'assets'), 'moleditpy_linux/assets'),
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
