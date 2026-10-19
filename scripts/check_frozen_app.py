"""Check the distributed executable version and dependency exclusions."""

import argparse
from pathlib import Path
import subprocess


def main() -> None:
    """Verify the actual executable and external Qt library layout."""
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    result = subprocess.run(
        [str(args.executable.resolve()), "--version"],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"MoleditPy {args.version}" in result.stdout, result.stdout + result.stderr
    contents = args.executable.parent
    if contents.name == "MacOS":
        contents = contents.parent
    bundled = [p for p in contents.rglob("*") if p.is_file()]
    forbidden = (
        "openbabel",
        "pyqt6",
        "qt6charts",
        "qtcharts",
        "qt6graphs",
        "qtgraphs",
        "qtdatavisualization",
        "qt6datavisualization",
        "qt6virtualkeyboard",
        "qtvirtualkeyboard",
    )
    assert not [str(p) for p in bundled if any(n in p.name.lower() for n in forbidden)]
    assert any("qt6core" in p.name.lower() or p.name == "QtCore" for p in bundled), (
        "Replaceable QtCore shared library missing"
    )
    print(result.stdout.strip())
    print(
        "Frozen application uses external Qt libraries and excludes retired dependencies."
    )


if __name__ == "__main__":
    main()
