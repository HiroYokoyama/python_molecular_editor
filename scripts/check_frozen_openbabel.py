"""Smoke-test a frozen MoleditPy executable and its isolated Open Babel child."""

import argparse
from pathlib import Path
import subprocess

from rdkit import Chem


def main() -> None:
    """Check version reporting and a real ethanol 3D conversion."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    executable = str(args.executable.resolve())
    version = subprocess.run(
        [executable, "--version"],
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )
    assert args.version in version.stdout, version.stdout + version.stderr
    molecule = Chem.MolFromSmiles("CCO")
    result = subprocess.run(
        [executable, "--openbabel-make3d"],
        input=Chem.MolToMolBlock(molecule),
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )
    converted = Chem.MolFromMolBlock(result.stdout, removeHs=False)
    assert converted is not None, result.stdout + result.stderr
    assert converted.GetNumAtoms() == 9
    assert converted.GetConformer().Is3D()
    print(f"MoleditPy {args.version}: frozen startup and Open Babel conversion passed")


if __name__ == "__main__":
    main()
