"""Smoke-test a frozen MoleditPy executable and its isolated Open Babel child."""

import argparse
import json
import math
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
    request = json.dumps(
        {
            "mol_block": result.stdout,
            "method": "UFF",
            "max_iters": 20,
            "chunk_size": 10,
        }
    )
    optimized = subprocess.run(
        [executable, "--openbabel-optimize"],
        input=request,
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )
    coordinates = json.loads(optimized.stdout)
    assert len(coordinates) == 9
    assert all(math.isfinite(value) for point in coordinates for value in point)
    print(f"MoleditPy {args.version}: frozen startup and Open Babel conversion passed")


if __name__ == "__main__":
    main()
