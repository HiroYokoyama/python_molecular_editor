# GitHub Actions Workflows

This directory contains the CI, test release, and production release workflows for MoleditPy.

## `tests.yml` - CI Tests

Runs on pushes and pull requests targeting `main` or `master`.

Test matrix:

| OS | Python versions |
| --- | --- |
| Ubuntu | 3.11, 3.12, 3.13, 3.14 |
| Windows | 3.14 |
| macOS | 3.14 |

The CI workflow installs the main package from `moleditpy/` and runs:

```bash
python tests/run_all_tests.py --no-cov --no-report --unit --integration
```

GUI tests are not run in CI.

## `test-release.yml` - Test Release

Manually builds a complete release candidate and publishes the Python packages to TestPyPI.

Run it from GitHub Actions:

1. Open **Actions**.
2. Select **Test Release**.
3. Click **Run workflow**.
4. Select the `main` branch.
5. Enter a test version such as `3.4.1.dev1`.

This workflow:

1. Validates that it is running on `main`.
2. Validates the version input.
3. Runs unit and integration tests on `windows-latest` with Python 3.12.
4. Updates `moleditpy/pyproject.toml` in the temporary runner checkout.
5. Runs `python scripts/sync_linux_version.py`.
6. Builds `MoleditPy` and `MoleditPy-linux` distributions.
7. Builds `dist/MoleditPy/MoleditPy.exe` with PyInstaller from the Linux package source.
8. Creates `MoleditPy_<version>_win64_portable.zip`.
9. Builds `MoleditPy_<version>_win64_setup.exe` with Inno Setup.
10. Uploads all built artifacts to the workflow run.
11. Publishes both Python packages to TestPyPI.

It does not commit, tag, push, or create a GitHub Release.

## `release.yml` - Production Release

Manually builds and publishes an official release.

Run it from GitHub Actions:

1. Open **Actions**.
2. Select **Release**.
3. Click **Run workflow**.
4. Select the `main` branch.
5. Enter a release version such as `3.4.1`.

This workflow:

1. Validates that it is running on `main`.
2. Validates the version input.
3. Rejects the release if the tag already exists.
4. Runs unit and integration tests on `windows-latest` with Python 3.12.
5. Updates `moleditpy/pyproject.toml`.
6. Runs `python scripts/sync_linux_version.py`.
7. Builds `MoleditPy` and `MoleditPy-linux` distributions.
8. Builds `dist/MoleditPy/MoleditPy.exe` with PyInstaller from the Linux package source.
9. Creates `MoleditPy_<version>_win64_portable.zip`.
10. Builds `MoleditPy_<version>_win64_setup.exe` with Inno Setup.
11. Commits the version and Linux sync changes.
12. Creates and pushes the release tag.
13. Publishes both Python packages to PyPI.
14. Creates a GitHub Release and attaches the distributions, installer, and portable ZIP.

## Required PyPI Setup

Both release workflows use PyPI trusted publishing with GitHub Actions OIDC.

Configure trusted publishers for these projects:

- `MoleditPy`
- `MoleditPy-linux`

For production PyPI, use workflow filename:

```text
release.yml
```

For TestPyPI, use workflow filename:

```text
test-release.yml
```

Recommended publisher settings:

| Field | Value |
| --- | --- |
| Owner | `HiroYokoyama` |
| Repository | `python_molecular_editor` |
| Workflow filename | `release.yml` or `test-release.yml` |
| Environment | empty, unless the workflow later adds one |

## Generated Artifacts

The release workflows generate:

- `dist/moleditpy/*.whl`
- `dist/moleditpy/*.tar.gz`
- `dist/moleditpy-linux/*.whl`
- `dist/moleditpy-linux/*.tar.gz`
- `dist/windows/MoleditPy_<version>_win64_setup.exe`
- `dist/windows/MoleditPy_<version>_win64_portable.zip`

## Files Modified During Release

`release.yml` commits these release-time changes:

- `moleditpy/pyproject.toml`
- `moleditpy-linux/`
- `windows-installer/script/script.iss`
- `windows-installer/script/MoleditPy.spec`

`test-release.yml` makes the same version and sync changes only inside the temporary GitHub Actions runner checkout.

## Windows Packaging Chain

Windows packaging uses the Linux package variant because Open Babel is disabled there.

Build order:

```text
moleditpy/pyproject.toml version update
-> scripts/sync_linux_version.py
-> moleditpy-linux package build
-> PyInstaller using windows-installer/script/MoleditPy.spec
-> dist/MoleditPy/MoleditPy.exe
-> portable ZIP from dist/MoleditPy/*
-> Inno Setup installer from dist/MoleditPy/*
```
