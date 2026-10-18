"""Entry point for Open Babel conversion inside a frozen application's child."""

import sys


def make3d_main() -> None:
    """Read one MOL block from stdin and write its generated 3D structure."""
    # Import through the worker to load RDKit first and initialize wheel data paths.
    from ..ui.calculation_worker import pybel

    if pybel is None:
        raise RuntimeError("Open Babel is not available.")
    molecule = pybel.readstring("mol", sys.stdin.read())
    molecule.addh()
    molecule.make3D()
    sys.stdout.write(molecule.write("mol"))
