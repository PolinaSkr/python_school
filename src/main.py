import distutils
import rdkit
import fastapi
from rdkit import Chem
from fastapi import FastAPI, HTTPException
import pandas as pd


def substructure_search(mols: list, mol: str):

    """
    :param mols: List of molecules as SMILES strings
    :param mol: Substructure as SMILES string
    :return: List of molecules as SMILES strings from mols that contain substructure mol
    """
    mol = Chem.MolFromSmiles(mol)
    return [m for m in mols if Chem.MolFromSmiles(m).HasSubstructMatch(mol)]


def add_molecule_to_db(identifier: str, smiles: str):
    if identifier in molecules_db:
        raise HTTPException(status_code=400, detail="Identifier already exists in database")
    molecules_db[identifier] = smiles
    return {"message": "Molecule successfully added in database"}

# res = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
# print(res)

app = FastAPI()

molecules_db = {}


@app.post("/add_molecule/", status_code=201,
          response_description="Add a new molecule in database", tags=["Database editing"])
def add_molecule(identifier: str, smiles: str):
    return add_molecule_to_db(identifier, smiles)


@app.get("/get_molecule/{identifier}", status_code=200,
         response_description="Get molecule by identifier", tags=["Database usage"])
def get_molecule(identifier: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found in database")
    return {"identifier": identifier, "smiles": molecules_db[identifier]}


@app.put("/update_molecule/{identifier}", status_code=200,
         response_description="Update molecule by identifier", tags=["Database editing"])
def update_molecule(identifier: str, smiles: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found in database")
    molecules_db[identifier] = smiles
    return {"message": "Molecule updated successfully"}


@app.delete("/delete_molecule/{identifier}", status_code=200,
            response_description="Delete molecule by identifier", tags=["Database editing"])
def delete_molecule(identifier: str):
    if identifier not in molecules_db:
        raise HTTPException(status_code=404, detail="Molecule not found in database")
    del molecules_db[identifier]
    return {"message": "Molecule deleted successfully"}


@app.get("/list_molecules/", status_code=200,
         response_description="List all molecules from database", tags=["Database usage"])
def list_molecules():
    if molecules_db:
        return molecules_db
    return {"message": "Database is empty"}


@app.post("/substructure_search/", status_code=200,
          response_description="Substructure search in database", tags=["Database usage"])
def substructure_search_api(mol: str):
    result = substructure_search(list(molecules_db.values()), mol)
    if result:
        return {"matching_molecules": result}
    return {"message": "No matching molecules in database"}


@app.post("/upload_molecules/", status_code=201,
          response_description="Upload molecules from file", tags=["Database editing"])
def upload_molecules(file_path: str):
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error reading file: {e}")

    identifier_list = df['identifier'].tolist()
    if len(identifier_list) != len(set(identifier_list)):
        return {"message": "Uploaded file contains duplicated identifiers. Uploading canceled."}

    smiles_list = df['smiles'].tolist()
    added_molecules = 0
    for ind, smiles in zip(identifier_list, smiles_list):
        try:
            add_molecule_to_db(ind, smiles)
            added_molecules += 1
        except HTTPException as e:
            if e.status_code == 400:
                continue
            else:
                raise e
    return {"message": f"{added_molecules} molecules from external database uploaded successfully"}
