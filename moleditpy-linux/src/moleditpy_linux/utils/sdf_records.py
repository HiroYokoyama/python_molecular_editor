#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MoleditPy — A Python-based molecular editing software

Author: Hiromichi Yokoyama
License: GPL-3.0 license
Repo: https://github.com/HiroYokoyama/python_molecular_editor
DOI: 10.5281/zenodo.17268532
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from io import BytesIO
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


@dataclass(frozen=True)
class SdfRecord:
    """One record of an SDF file; ``mol`` is None when RDKit cannot read it."""

    index: int
    mol: Optional[Chem.Mol]
    name: str

    @property
    def readable(self) -> bool:
        """Whether RDKit could build a molecule from this record."""
        return self.mol is not None

    @property
    def formula(self) -> str:
        """Molecular formula, or an empty string for an unreadable record."""
        if self.mol is None:
            return ""
        try:
            return str(CalcMolFormula(self.mol))
        except (RuntimeError, ValueError):
            return ""

    @property
    def num_atoms(self) -> int:
        """Atom count including explicit hydrogens, 0 when unreadable."""
        return 0 if self.mol is None else int(self.mol.GetNumAtoms())

    def display_name(self) -> str:
        """The record's title line, or a numbered placeholder when it is blank."""
        return self.name or f"Molecule {self.index + 1}"


def read_sdf_records(sdf_text: str) -> List[SdfRecord]:
    """Parse every record of already-decoded SDF text, in file order.

    The text is fed to RDKit as a stream rather than a path, so the caller's
    encoding fallbacks apply and each record keeps its data fields. A record
    RDKit rejects stays in the list with ``mol`` None, so indices match the
    file and the user can see that it was skipped.
    """
    # Trailing blank lines after the last "$$$$" would read as one more
    # (empty, unreadable) record.
    text = sdf_text.rstrip() + "\n"
    titles = _record_titles(text)
    records: List[SdfRecord] = []
    supplier = Chem.ForwardSDMolSupplier(BytesIO(text.encode("utf-8")), removeHs=False)
    index = 0
    while True:
        try:
            mol = next(supplier)
        except StopIteration:
            break
        except (RuntimeError, ValueError) as e:
            logging.warning("SDF record %d could not be read: %s", index + 1, e)
            mol = None
        if mol is not None:
            name = mol.GetProp("_Name").strip()
        else:
            # RDKit gives nothing back for a record it rejects; take the title
            # line from the text so the list still names it.
            name = titles[index] if index < len(titles) else ""
        records.append(SdfRecord(index=index, mol=mol, name=name))
        index += 1
    return records


def _record_titles(sdf_text: str) -> List[str]:
    """First line (the title) of each "$$$$"-separated record, stripped."""
    titles: List[str] = []
    at_start = True
    for line in sdf_text.splitlines():
        if at_start:
            titles.append(line.strip())
            at_start = False
        if line.strip() == "$$$$":
            at_start = True
    return titles
