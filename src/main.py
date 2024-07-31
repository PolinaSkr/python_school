import rdkit
from rdkit import Chem


def substructure_search(mols: list, mol: str):

    """
    :param mols: List of molecules as SMILES strings
    :param mol: Substructure as SMILES string
    :return: List of molecules as SMILES strings from mols that contain substructure mol
    """
    mol = Chem.MolFromSmiles(mol)
    return [m for m in mols if Chem.MolFromSmiles(m).HasSubstructMatch(mol)]


res = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
print(res)
