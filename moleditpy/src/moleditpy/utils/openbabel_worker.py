"""Native Open Babel operations, loaded only in isolated child processes."""

import json
import os
from pathlib import Path
import sys
from typing import Any

# Load RDKit before Open Babel to avoid native symbol conflicts on Linux.
from rdkit import Chem


def _load_pybel() -> Any:
    """Load Open Babel and find the parameter files bundled with its wheel."""
    Chem.GetPeriodicTable()
    import openbabel
    from openbabel import pybel

    datadir = os.environ.get("BABEL_DATADIR")
    if not datadir or not (Path(datadir) / "UFF.prm").is_file():
        ob_base = Path(openbabel.__file__).parent
        candidates = [ob_base / "bin" / "data"]
        candidates.extend((ob_base / "share" / "openbabel").glob("*"))
        for candidate in candidates:
            if (candidate / "UFF.prm").is_file():
                os.environ["BABEL_DATADIR"] = str(candidate)
                break
    return pybel


def make3d_main() -> None:
    """Read one MOL block from stdin and write its generated 3D structure."""
    pybel = _load_pybel()
    molecule = pybel.readstring("mol", sys.stdin.read())
    molecule.addh()
    molecule.make3D()
    sys.stdout.write(molecule.write("mol"))


def optimize_main() -> None:
    """Read an optimization request and return coordinates in atom order."""
    request = json.load(sys.stdin)
    pybel = _load_pybel()
    molecule = pybel.readstring("mol", request["mol_block"])
    forcefield = pybel.ob.OBForceField.FindForceField(request["method"].lower())
    if not forcefield:
        raise ValueError(f"Forcefield '{request['method']}' not found.")
    if not forcefield.Setup(molecule.OBMol):
        raise RuntimeError(f"Failed to load parameters for '{request['method']}'.")
    forcefield.SteepestDescent(100)
    for _ in range(0, request["max_iters"], request["chunk_size"]):
        forcefield.ConjugateGradients(request["chunk_size"])
    forcefield.GetCoordinates(molecule.OBMol)
    coordinates = []
    for i in range(1, molecule.OBMol.NumAtoms() + 1):
        atom = molecule.OBMol.GetAtom(i)
        coordinates.append([atom.GetX(), atom.GetY(), atom.GetZ()])
    json.dump(coordinates, sys.stdout)
