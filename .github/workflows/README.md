# GitHub Actions workflows

The YAML files are the source of truth for triggers, matrices, and publishing.
All desktop platforms use the single `moleditpy/` package. Python 3.9–3.14 is
supported, with PySide6-Essentials (Qt 6.9 on Python 3.9; Qt 6.10–6.11 otherwise).
Open Babel and PyQt6 are excluded from the core and frozen bundles.

| Workflow | Purpose |
| --- | --- |
| `tests.yml` | Unit, integration, E2E and GUI tests, plus pinned lint/type checks. Linux tests Python 3.9–3.14; Windows and macOS test 3.9, 3.13 and 3.14. Windows/Python 3.13 measures coverage. |
| `full-gui-tests.yml` | Native VTK/PyVista rendering suites on Linux/Python 3.11 and 3.13 and Windows/Python 3.13. Launch probes cover all three OSes and Python 3.9–3.14. Linux and macOS require a working probe; Windows launch and rendering tests can skip when OpenGL is unavailable. macOS runs the isolated subprocess launch test. |
| `windows-exe-smoke.yml` | Builds the actual main-package Windows executable and checks its version, replaceable Qt libraries, and excluded dependencies. Uploads the tested executable; does not publish a package. |
| `test-pip-install.yml` | Manual PyPI/TestPyPI install validation across all three OSes and Python 3.9–3.14 (18 combinations). Specify `5.0.0a1` to test the alpha; leaving the version empty selects a stable release. Runs E2E tests with `MOLEDITPY_USE_INSTALLED=1`. |
| `macos-app.yml` | Reusable/manual Apple Silicon app build, tests, frozen-app validation, ad-hoc signing, and ZIP artifact creation. Can attach an artifact to an existing release when explicitly requested. |
| `test-release.yml` | Manual packaging rehearsal from `main`, with Windows/macOS artifacts and publication to TestPyPI. Does not create the production release or push a production tag. Use a new TestPyPI version for each upload. |
| `release.yml` | Manual production PyPI and GitHub publication; stable releases run from `main`, prereleases may run from another branch. |
| `zenodo.yml` | Separate manual Zenodo publication workflow. |
| `test-zenodo.yml` | Manual Zenodo Sandbox workflow. The `draft` input controls whether the sandbox record is published. |

## Manual plugin tests

Plugin test workflows that install MoleditPy expose a manual `moleditpy_ref`
input for testing against an app branch, tag, or full commit SHA. Leaving it blank
preserves normal CI selection: the repository's `MOLEDITPY_REF` variable, when configured, otherwise
the migration branch `migrate/pyside6-without-openbabel`. Select the plugin's own
branch separately in the workflow dispatch branch selector. A published app tag
such as `v5.0.0a1` selects that release's source for plugin compatibility tests.

The main app's `test-pip-install.yml` already accepts an exact published version
and PyPI/TestPyPI source. This is separate from checking out an app source ref.

## Running tests locally

```sh
python -m pip install -e moleditpy/
python -m pip install pytest pytest-qt pytest-mock pytest-cov pytest-timeout
python tests/run_all_tests.py --pylint
```

Use `QT_QPA_PLATFORM=offscreen`, `QT_API=pyside6`, and
`PYTEST_QT_API=pyside6` for the regular headless suites. The no-flag runner
includes unit, integration, E2E and GUI tests. Native full-GUI tests require
`--full-gui` and a working display/OpenGL context; Linux CI uses xvfb/Mesa.
`MOLEDITPY_USE_INSTALLED=1` selects the installed package for E2E and full-GUI
validation instead of prepending the checkout's source directory.

## Production release

1. Set the version in `moleditpy/pyproject.toml`, update the changelog, and pass CI.
2. Dispatch `release.yml` from the intended branch with that exact version.
3. Stable versions require `main`. Another branch must use an alpha/beta/rc
   version such as `5.0.0a1`. The tag must not already exist.
4. Choose the Zenodo input deliberately. Alpha/beta/rc releases skip Zenodo
   regardless of that input and are not marked as GitHub's latest stable release.

The workflow validates the version/ref, builds and smoke-tests the macOS app,
then runs Windows tests and builds the wheel/sdist, executable, portable ZIP,
and installer. Frozen checks validate version reporting, external Qt libraries,
and exclusion of Open Babel, PyQt6, and unused Qt add-ons; they do not perform
an Open Babel conversion. Native rendering is checked by `full-gui-tests.yml`.

After both platform builds succeed, the workflow pushes a new tag, publishes the
main package to PyPI from a Linux job using trusted publishing, and creates the
GitHub release. A prerelease tags its dispatched commit without a version-bump
commit or branch push. A stable release may commit updated version/download
metadata. The optional final Zenodo job runs only for stable releases.

The shortcut helper and plugin repositories have separate publication processes.
Releasing the main app does not publish their migration branches or alter the
retired Linux project's existing PyPI downloads.
