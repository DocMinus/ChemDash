#!/usr/bin/env python
# -*- coding: utf-8 -*-
# changed do work without files, now SVG only in memory
# -22 red; -233 green

from rdkit import Chem
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.CRITICAL)


from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict


def elemental_composition(smiles: str) -> dict:
    """Get Atomic counts, including hydrogen atoms, and any charge.
    Taken from Rdkit discussions on Github
    :param molecule: The molecule to analyze
    :rtype: A dictionary.
    """

    comp = defaultdict(lambda: 0)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol_ok = True

    try:
        Chem.SanitizeMol(mol)
    except ValueError as _e:
        mol_ok = False
        pass

    if mol_ok and len(smiles) > 0:
        mol = Chem.AddHs(mol)

        # Get atom counts
        for atom in mol.GetAtoms():
            comp[atom.GetAtomicNum()] += 1

        # If charged, add charge as "atomic number" 0
        charge = Chem.GetFormalCharge(mol)
        if charge != 0:
            comp[0] = charge

    return comp


def oxy_balance(smiles: str) -> str:

    composition = elemental_composition(smiles)
    o_balance = 0
    if composition:
        o_count = composition[8]
        c_count = composition[6]
        h_count = composition[1]
        m_count = 0  # metal count; most cases are 0. but could be calculated.
        weight = Descriptors.MolWt(Chem.MolFromSmiles(smiles, sanitize=True))
        o_balance = (-1600 / (weight)) * (
            (2 * c_count) + (h_count / 2) + m_count - o_count
        )

    return f"{round(o_balance,2)}"


def n_count(smiles: str) -> int:
    composition = elemental_composition(smiles)
    n_count = 0
    if composition:
        n_count = composition[7]
    return n_count


# print(oxy_balance("C1(C=C([N+](=O)[O-])C(C)=C([N+]([O-])=O)C=1)[N+]([O-])=O"))
# example TNT which should give -73.9
