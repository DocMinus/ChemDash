# Updated: Corrected spelling of class name & renamed file for clarity.
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.logger().setLevel(RDLogger.CRITICAL)

from .chemtools import RDKIT_DESCRIPS, RDKIT_DESCRIPS_HEADERS, clean_smiles


def get_clean_ids(_mollist: list) -> list:
    """
    returns a list of molecule IDs used at the end of the property calculations.
    assumes that an ID exists. #TODO error checking, low prio. (simple if then should do).
    See comment in rdkit_calc: not necessary for a single mol, but prepares for future multiple mols at once.

    :param _list: molecule list (after cleaning, but could be used also before)
    :return: mol IDs as list
    """
    _idlist = []
    for mol in _mollist:
        name = mol.GetProp("_Name")
        _idlist.append(name)
    return _idlist


def rdkit_calc(cleaned_mol_list: list) -> pd.DataFrame:
    """
    calculates RDkit descriptors. Based on a list.
    Not necessary for the single mole cases in this proof of concept,
    but can be easily expanded if and when multiple mols at once are to be considered.

    :param cleaned_mol_list: list of RDkit molecule objects
    :return: pandas df containing calculated properties (no structures, but now with ID)
    """

    calcrdk = [RDKIT_DESCRIPS.CalcDescriptors(mol) for mol in cleaned_mol_list]
    transposed_properties = zip(*calcrdk)
    headers = RDKIT_DESCRIPS_HEADERS
    named_properties = dict(zip(headers, transposed_properties))
    properties_df1 = pd.DataFrame(data=named_properties)

    properties_ids_df = pd.DataFrame(
        data=get_clean_ids(cleaned_mol_list), columns=["cleanIDs"]
    )
    return pd.concat([properties_ids_df, properties_df1], axis=1)


class Molecule:
    """
    Simple clean and convert of smiles to rdkit object
    Adding more stuff possible, but not necessary for this proof of concept.

    :param molecule: mol as smiles \n
    :rtype: cleanmolobj() = rdkit object \n
            cleansmiles() = rkdit object as smiles \n
            sumformula() \n
            molwt() \n
            moldescriptors() = rdkit calculated descriptors
    """

    def __init__(self, schmiles: str):
        self._schmiles = schmiles
        self._smi_df = pd.DataFrame({"ID": ["dummy"], "molecule": [self._schmiles]})
        self._in_mols_df = self._smi_df[
            list(self._smi_df.columns[1:2]) + list(self._smi_df.columns[0:1])
        ]
        self._mol_list = clean_smiles(self._in_mols_df)
        self._mol = self._mol_list[0]

    def clean_molobj(self):
        return self._mol

    def cleaned_smiles(self):
        return Chem.rdmolfiles.MolToSmiles(self._mol)

    def sumformula(self):
        return Chem.rdMolDescriptors.CalcMolFormula(self._mol)

    def molwt(self):
        return Descriptors.MolWt(self._mol)

    def moldescriptors(self):
        return rdkit_calc([self._mol])
