"""Validate project contents before replacing the live document."""

from __future__ import annotations

import base64
import io
import math
import pickle
import pickletools
from typing import Any, BinaryIO

import numpy as np
from rdkit import Chem

from .molecular_data import PointTuple


def _integer(value: Any, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError("Project integer field is invalid")
    return value


def _coordinate(value: Any) -> None:
    try:
        if isinstance(value, bool) or not math.isfinite(float(value)):
            raise ValueError("Project coordinates must be finite numbers")
    except OverflowError as exc:
        raise ValueError("Project coordinate is out of range") from exc


def _validate_2d(structure: Any) -> None:
    if not isinstance(structure, dict):
        raise ValueError("2D structure must be an object")
    atoms, bonds = structure.get("atoms", []), structure.get("bonds", [])
    if not isinstance(atoms, list) or not isinstance(bonds, list):
        raise ValueError("Project atoms and bonds must be lists")
    ids: set[int] = set()
    for atom in atoms:
        atom_id = _integer(atom["id"])
        if atom_id in ids or not isinstance(atom["symbol"], str) or not atom["symbol"]:
            raise ValueError("Project atoms need unique IDs and element symbols")
        ids.add(atom_id)
        _coordinate(atom["x"])
        _coordinate(atom["y"])
        _integer(atom.get("charge", 0), -(2**31))
        if atom.get("charge", 0) >= 2**31:
            raise ValueError("Project charge is out of range")
        _integer(atom.get("radical", 0))
        if atom.get("rs") not in (None, "R", "S"):
            raise ValueError("Project R/S label must be R or S")
    seen: set[frozenset[int]] = set()
    for bond in bonds:
        first, second = _integer(bond["atom1"]), _integer(bond["atom2"])
        pair = frozenset((first, second))
        if first == second or first not in ids or second not in ids or pair in seen:
            raise ValueError("Project bond has invalid atom references")
        seen.add(pair)
        if isinstance(bond["order"], bool) or bond["order"] not in (1, 1.5, 2, 3):
            raise ValueError("Project bond order is invalid")
        if _integer(bond.get("stereo", 0)) not in (0, 1, 2, 3, 4):
            raise ValueError("Project bond stereochemistry is invalid")
    minimum = max(ids) + 1 if ids else 0
    _integer(structure.get("next_atom_id", minimum), minimum)


def _validate_molecule(binary: Any) -> Chem.Mol:
    mol = Chem.Mol(binary)
    for conformer in mol.GetConformers():
        if not np.isfinite(conformer.GetPositions()).all():
            raise ValueError("Project 3D coordinates must be finite")
    return mol


def validate_project_data(data: dict[str, Any]) -> None:
    """Reject corrupt geometry without changing application or plugin state."""
    if "2d_structure" in data:
        _validate_2d(data["2d_structure"])
    structure = data.get("3d_structure")
    if structure is None:
        return
    if not isinstance(structure, dict):
        raise ValueError("3D structure must be an object")
    atoms = structure.get("atoms", [])
    if not isinstance(atoms, list):
        raise ValueError("3D atoms must be a list")
    encoded = structure.get("mol_binary_base64")
    if encoded:
        binary: Any = base64.b64decode(encoded.encode("ascii"), validate=True)
        mol = _validate_molecule(binary)
        if atoms and len(atoms) != mol.GetNumAtoms():
            raise ValueError("3D atom metadata does not match its molecule")
    seen: set[int] = set()
    for atom in atoms:
        index = _integer(atom["index"])
        if index >= len(atoms) or index in seen:
            raise ValueError("3D atom index is invalid")
        seen.add(index)
        for axis in ("x", "y", "z"):
            _coordinate(atom.get(axis, 0.0))


# The only non-builtin type a raw snapshot contains: 2D atom positions are
# PointTuple, a plain tuple subclass whose reconstruction runs no user code.
_RAW_PROJECT_CLASSES: dict[tuple[str, str], type] = {
    ("moleditpy.core.molecular_data", "PointTuple"): PointTuple,
}


class _RestrictedUnpickler(pickle.Unpickler):
    """Unpickler that refuses every global reference except PointTuple.

    A raw project snapshot is otherwise built from dicts, lists, tuples,
    strings, bytes and numbers, none of which needs ``find_class``. Any other
    reference is a pickle trying to import a callable, which is the
    arbitrary-code-execution vector of the format.
    """

    def find_class(self, module: str, name: str) -> Any:
        try:
            return _RAW_PROJECT_CLASSES[(module, name)]
        except KeyError:
            raise pickle.UnpicklingError(
                f"Raw project files may not reference {module}.{name}"
            ) from None


# Everything a snapshot needs: data, containers, memo and framing, plus the
# global lookup and NEWOBJ that rebuild a PointTuple. Left out on purpose are
# the opcodes that call or configure an object (REDUCE, BUILD, INST, OBJ,
# NEWOBJ_EX), persistent IDs, out-of-band buffers and the EXT opcodes: a cached
# copyreg extension is returned without ever reaching find_class.
_RAW_PROJECT_OPCODES = frozenset(
    {
        "PROTO", "FRAME", "STOP", "MARK", "POP", "POP_MARK", "DUP",
        "MEMOIZE", "GET", "BINGET", "LONG_BINGET", "PUT", "BINPUT", "LONG_BINPUT",
        "NONE", "NEWTRUE", "NEWFALSE",
        "INT", "BININT", "BININT1", "BININT2", "LONG", "LONG1", "LONG4",
        "FLOAT", "BINFLOAT",
        "STRING", "BINSTRING", "SHORT_BINSTRING",
        "UNICODE", "SHORT_BINUNICODE", "BINUNICODE", "BINUNICODE8",
        "BINBYTES", "SHORT_BINBYTES", "BINBYTES8", "BYTEARRAY8",
        "EMPTY_LIST", "LIST", "APPEND", "APPENDS",
        "EMPTY_TUPLE", "TUPLE", "TUPLE1", "TUPLE2", "TUPLE3",
        "EMPTY_DICT", "DICT", "SETITEM", "SETITEMS",
        "EMPTY_SET", "ADDITEMS", "FROZENSET",
        "GLOBAL", "STACK_GLOBAL", "NEWOBJ",
    }
)  # fmt: skip


def load_raw_project(stream: BinaryIO) -> Any:
    """Unpickle a .pmeraw stream without allowing any imports or calls."""
    data = stream.read()
    try:
        for opcode, _argument, _position in pickletools.genops(data):
            if opcode.name not in _RAW_PROJECT_OPCODES:
                raise pickle.UnpicklingError(
                    f"Raw project files may not use the {opcode.name} pickle opcode"
                )
    except (ValueError, OverflowError, IndexError, TypeError) as exc:
        raise pickle.UnpicklingError(f"Malformed raw project file: {exc}") from exc
    return _RestrictedUnpickler(io.BytesIO(data)).load()


def validate_raw_project_data(data: Any) -> None:
    """Validate legacy pickle geometry before any live state is cleared."""
    if not isinstance(data, dict):
        raise ValueError("Raw project must be an object")
    atoms, bonds = data.get("atoms", {}), data.get("bonds", {})
    if not isinstance(atoms, dict) or not isinstance(bonds, dict):
        raise ValueError("Raw project atoms and bonds must be mappings")
    structure = {
        "atoms": [
            dict(atom, id=atom_id, x=atom["pos"][0], y=atom["pos"][1])
            for atom_id, atom in atoms.items()
        ],
        "bonds": [
            dict(bond, atom1=key[0], atom2=key[1]) for key, bond in bonds.items()
        ],
        "next_atom_id": data.get("_next_atom_id", max(atoms) + 1 if atoms else 0),
    }
    _validate_2d(structure)
    if data.get("mol_3d") is not None:
        _validate_molecule(data["mol_3d"])
