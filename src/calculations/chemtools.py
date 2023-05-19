"""
Reduced chemtools module for use in ChemDash proof of concept based on
Version 0.5.6 (May 24, 07:30:00 2021)
Update: 2023-01-02
@author: Alexander Minidis (DocMinus)
Copyright (c) 2021-2023 DocMinus
"""

import re

import pandas as pd

# RDkit stuff
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

# Important, or else waaaay too many RDkit details in output
RDLogger.logger().setLevel(RDLogger.CRITICAL)


def is_smiles(text: str) -> bool:
    """
    Determine if input possibly is a smiles or just a regular text
    Tokenize the input and compare the input length vs token length
    If not the same, then a regular word or phrase.
    based on: https://github.com/pschwllr/MolecularTransformer
        Input:
            string
        Output:
            boolean, True if smiles, False if regular text
    """
    # This simple pre-check seems to work for catching plain numbers as ID (else detected as smiles)
    # normally isnumeric could be used, but if it is a number, it's from numpy and throws error.
    if not isinstance(text, str):
        return False

    pattern = "(\[[^\]]+]|Si|Ti|Al|Zn|Pd|Pt|Cu|Br?|Cl?|N|O|S|P|F|I|B|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(text)]
    length = len("".join(tokens))
    # empty could give len 0!
    # >= probably not necessary, but safer than ==
    if length and length >= len(text):
        return True
    else:
        return False


# deleted stuff

""" ######## Cleaning of smiles ########
 3 (well, 4) functions in total:
 * clean_smiles - pandas of molecules stemming from input smiles. The minimum cleaner.
 * clean_all: wrapper combining clean_smiles and step2, 3, 4. Use for thoroughness or nasty datasources
 * wrapper allows for flexibility in only using one func at a time
"""


def clean_smiles(mols_raw: pd.DataFrame) -> list:
    """
    First round of cleaning of the input smiles structures.
    Creates RDKit object.
    Cleaning includes some standard normalization.
    Removes completely faulty sh** right from the start, works well for relatively good sources

    :param mols_raw: either dictionary (deprecated, though still in code), or pandas, 'preformatted' order
    :return: list containing rdkit molobjects
    """
    # NOTE: check here for details: https://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization

    if not isinstance(mols_raw, dict) and not isinstance(mols_raw, pd.DataFrame):
        raise TypeError("Dictionary or pandas DataFrame expected")

    if isinstance(mols_raw, pd.DataFrame):
        mols = pd.Series(
            mols_raw.iloc[:, 0].values, index=mols_raw.iloc[:, 1]
        ).to_dict()
    else:
        mols = mols_raw

    num_entries = len(mols)
    molecule_list = []
    for name, smi in mols.items():
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            continue
        try:
            Chem.SanitizeMol(mol)

        except ValueError as _e:
            print("ID: ", name, " ", _e)
            continue

        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=(
                Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES
            ),
        )
        mol = rdMolStandardize.Normalize(mol)
        # N.B. it seems that in first iteration strings as digits still are seen as
        # integers which throws an error from rdkit. this an additional conversion to str.
        # more readable, though I wonder if one could just write enforced name=str(name)?
        if not isinstance(name, str):
            name = str(name)
        mol.SetProp("_Name", name)
        molecule_list.append(mol)

    if len(molecule_list) == 0:  # "if not molecule_list:" should also work
        raise ImportError("!!! No molecules to work with!!!")
    elif len(molecule_list) < num_entries:
        print(
            "!N.B.: ",
            num_entries - len(molecule_list),
            " faulty entries have been filtered out!",
        )

    return molecule_list


def clean_step2(mol_object: list) -> list:
    """
    Neutralize, Reionize, Keep parent fragment & clean some metal bonds
    of the rdkit object

    :param mol_object: list of rdkit objects
    return: new list of rdkit objects
    """
    uncharger = rdMolStandardize.Uncharger()
    disconnector = rdMolStandardize.MetalDisconnector()
    molecule_list = []
    for mol in mol_object:
        name = mol.GetProp("_Name")
        mol = uncharger.uncharge(mol)

        try:
            uncharger.uncharge(rdMolStandardize.FragmentParent(mol))

        except ValueError as _e:
            print("ID: ", name, " ", _e)
            continue

        mol = uncharger.uncharge(rdMolStandardize.FragmentParent(mol))
        mol = disconnector.Disconnect(mol)
        # Fragmentor seems to lose the name thus need to save and set manually
        mol.SetProp("_Name", name)
        molecule_list.append(mol)

    if len(molecule_list) < len(mol_object):
        print(
            "!N.B.: ",
            len(mol_object) - len(molecule_list),
            " faulty entries have been filtered out!",
        )

    return molecule_list


def clean_step3(mol_object: list) -> list:
    """
    Neutralizes charges in a rdkit-object (containing list)
    """
    uncharger = rdMolStandardize.Uncharger()
    molecule_list = []
    for mol in mol_object:
        mol = uncharger.uncharge(mol)
        molecule_list.append(mol)

    return molecule_list


def clean_step4(mol_object: list) -> list:
    """
    Reionizes and keeps parent fragment of rdkit object(containing list)
    """
    uncharger = rdMolStandardize.Uncharger()
    md = rdMolStandardize.MetalDisconnector()
    molecule_list = []
    for mol in mol_object:
        # name = mol.GetProp("_Name")
        # using name here seems to mess things up;
        # on the other hand, not necessary anyway since rdkit takes care of it at this stage
        mol = uncharger.uncharge(rdMolStandardize.FragmentParent(mol))
        mol = md.Disconnect(mol)

        # mol.SetProp("_Name", name)
        molecule_list.append(mol)

    return molecule_list


def clean_all(mols: pd.DataFrame) -> list:
    """
    Combines all cleaning functions into one. Could four steps cleaning be done in one step?
    Yes, but this is more flexible, can in principle skip one or more steps.
    Use for thoroughness or for nasty datasources

    :param mols: pandas df containing smiles strings
    :return: list of rdkit objects
    """

    _x0 = clean_smiles(mols)
    _x1 = clean_step2(_x0)
    _x2 = clean_step3(_x1)
    _x3 = clean_step4(_x2)

    return _x3


# deleted stuff

# VARIABLES
RDKIT_DESCRIPS = MolecularDescriptorCalculator(
    [
        "MolLogP",
        "MolMR",
        "MolWt",
        "NHOHCount",
        "NOCount",
        "FractionCSP3",
        "RingCount",
        "NumAliphaticCarbocycles",
        "NumAliphaticHeterocycles",
        "NumAliphaticRings",
        "NumAromaticCarbocycles",
        "NumAromaticHeterocycles",
        "NumAromaticRings",
        "NumHAcceptors",
        "NumHDonors",
        "NumHeteroatoms",
        "NumRadicalElectrons",
        "NumRotatableBonds",
        "NumSaturatedCarbocycles",
        "NumSaturatedHeterocycles",
        "NumSaturatedRings",
        "NumValenceElectrons",
        "HeavyAtomCount",
        "HeavyAtomMolWt",
        "TPSA",
        "MaxEStateIndex",
        "BalabanJ",
        "BertzCT",
        "Chi0",
        "Chi0n",
        "Chi0v",
        "Chi1",
        "Chi1n",
        "Chi1v",
        "Chi2n",
        "Chi2v",
        "Chi3n",
        "Chi3v",
        "Chi4n",
        "Chi4v",
        "EState_VSA1",
        "EState_VSA10",
        "EState_VSA11",
        "EState_VSA2",
        "EState_VSA3",
        "EState_VSA4",
        "EState_VSA5",
        "EState_VSA6",
        "EState_VSA7",
        "EState_VSA8",
        "EState_VSA9",
        "ExactMolWt",
        "FpDensityMorgan1",
        "FpDensityMorgan2",
        "FpDensityMorgan3",
        "HallKierAlpha",
        "Ipc",
        "Kappa1",
        "Kappa2",
        "Kappa3",
        "LabuteASA",
        "MaxAbsEStateIndex",
        "MaxAbsPartialCharge",
        "MaxPartialCharge",
        "MinAbsEStateIndex",
        "MinAbsPartialCharge",
        "MinEStateIndex",
        "MinPartialCharge",
        "PEOE_VSA1",
        "PEOE_VSA10",
        "PEOE_VSA11",
        "PEOE_VSA12",
        "PEOE_VSA13",
        "PEOE_VSA14",
        "PEOE_VSA2",
        "PEOE_VSA3",
        "PEOE_VSA4",
        "PEOE_VSA5",
        "PEOE_VSA6",
        "PEOE_VSA7",
        "PEOE_VSA8",
        "PEOE_VSA9",
        "qed",
        "SlogP_VSA1",
        "SlogP_VSA10",
        "SlogP_VSA11",
        "SlogP_VSA12",
        "SlogP_VSA2",
        "SlogP_VSA3",
        "SlogP_VSA4",
        "SlogP_VSA5",
        "SlogP_VSA6",
        "SlogP_VSA7",
        "SlogP_VSA8",
        "SlogP_VSA9",
        "SMR_VSA1",
        "SMR_VSA10",
        "SMR_VSA2",
        "SMR_VSA3",
        "SMR_VSA4",
        "SMR_VSA5",
        "SMR_VSA6",
        "SMR_VSA7",
        "SMR_VSA8",
        "SMR_VSA9",
        "VSA_EState1",
        "VSA_EState10",
        "VSA_EState2",
        "VSA_EState3",
        "VSA_EState4",
        "VSA_EState5",
        "VSA_EState6",
        "VSA_EState7",
        "VSA_EState8",
        "VSA_EState9",
    ]
)

RDKIT_DESCRIPS_HEADERS = list(RDKIT_DESCRIPS.GetDescriptorNames())
