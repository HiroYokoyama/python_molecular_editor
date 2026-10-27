"""Conservative formal-charge repair for drawn molecular graphs."""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, Set

from rdkit import Chem, rdBase

if TYPE_CHECKING:
    from .molecular_data import MolecularData


def infer_charge_adjustments(data: MolecularData) -> Dict[int, int]:
    """Return unambiguous repairs without changing the editor's graph.

    Keep RDKit-valid charges. For invalid B/C/N/O valences, accept a repair
    only when exactly one of -1, 0, +1 passes RDKit's strict valence check.
    Missing hydrogens do not establish a charge: a singly bonded neutral O,
    for example, may be an alcohol and is left alone. Radicals, other elements,
    and atoms adjacent to unsupported symbols are outside this heuristic.
    """
    mol = Chem.RWMol()
    indices = {}
    bond_types = {
        1.0: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
    }
    with rdBase.BlockLogs():
        for atom_id, record in data.atoms.items():
            try:
                atom = Chem.Atom(record["symbol"])
            except (RuntimeError, ValueError):
                atom = Chem.Atom(0)
            atom.SetFormalCharge(record.get("charge", 0))
            atom.SetNumRadicalElectrons(record.get("radical", 0))
            indices[atom_id] = mol.AddAtom(atom)

        excluded: Set[int] = set()
        for (first, second), record in data.bonds.items():
            bond_type = bond_types.get(float(record["order"]))
            if first not in indices or second not in indices or bond_type is None:
                excluded.update((first, second))
                continue
            mol.AddBond(indices[first], indices[second], bond_type)

        changes = {}
        for atom_id, idx in indices.items():
            atom = mol.GetAtomWithIdx(idx)
            if (
                atom_id in excluded
                or atom.GetAtomicNum() not in (5, 6, 7, 8)
                or atom.GetNumRadicalElectrons()
                or any(n.GetAtomicNum() == 0 for n in atom.GetNeighbors())
            ):
                continue
            try:
                atom.UpdatePropertyCache(strict=True)
                continue
            except (RuntimeError, ValueError):
                pass

            candidates = []
            for charge in (-1, 0, 1):
                atom.SetFormalCharge(charge)
                try:
                    atom.UpdatePropertyCache(strict=True)
                    candidates.append(charge)
                except (RuntimeError, ValueError):
                    pass
            atom.SetFormalCharge(data.atoms[atom_id].get("charge", 0))
            if len(candidates) == 1:
                changes[atom_id] = candidates[0]
    return changes
