# -*- mode: python ; coding: utf-8 -*-

import os
from pathlib import Path


spec_dir = Path(SPECPATH).resolve()
repo_root = spec_dir.parent.parent
linux_package_dir = repo_root / 'moleditpy-linux' / 'src' / 'moleditpy_linux'
src_dir = linux_package_dir.parent  # moleditpy-linux/src/ — makes `import moleditpy_linux` work

app_version = os.environ.get('MOLEDITPY_APP_VERSION', '0.0.0')

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
        # Frozen apps carry no dist metadata, so constants._get_version() falls
        # back to walking up for a pyproject.toml — ship one where it looks.
        (str(linux_package_dir.parent.parent / 'pyproject.toml'), 'moleditpy_linux'),
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
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    # Left off deliberately: AppleEvent emulation can stall a bundle launched
    # from a shell, and the bundle declares no document types anyway.
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='MoleditPy',
)
app = BUNDLE(
    coll,
    name='MoleditPy.app',
    icon=str(linux_package_dir / 'assets' / 'icon.icns'),
    bundle_identifier='io.github.hiroyokoyama.moleditpy',
    version=app_version,
    info_plist={
        'CFBundleShortVersionString': app_version,
        'CFBundleVersion': app_version,
        'NSHighResolutionCapable': True,
        'LSMinimumSystemVersion': '11.0',
        'NSHumanReadableCopyright': 'Copyright (C) Hiromichi Yokoyama. GPL-3.0.',
    },
)
