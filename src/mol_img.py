#!/usr/bin/env python
# -*- coding: utf-8 -*-
# changed do work without files, now SVG only in memory

import base64
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.CRITICAL)


def mol_image(smiles: str) -> base64:
    """returns a base64 based SVG image from smiles string, in memory, no tmp file"""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol_ok = True
    _image = ""

    try:
        Chem.SanitizeMol(mol)
    except ValueError as _e:
        mol_ok = False
        pass

    if mol_ok and len(smiles) > 0:
        rdDepictor.SetPreferCoordGen(True)
        _d = rdMolDraw2D.MolDraw2DSVG(200, 200)
        rdMolDraw2D.SetDarkMode(_d)
        # _d.drawOptions().useBWAtomPalette()
        _d.drawOptions().padding = 0.1
        _d.drawOptions().scalingFactor = 100
        rdMolDraw2D.PrepareAndDrawMolecule(_d, mol)
        _d.FinishDrawing()
        _image = _d.GetDrawingText()

    encoded_svgimg = _image.encode("ascii")
    return base64.b64encode(encoded_svgimg).decode()
