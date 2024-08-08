import pytest
from fastapi.testclient import TestClient
from main import app, substructure_search, add_molecule_to_db, molecules_db

client = TestClient(app)

# Clear the molecules_db before each test
@pytest.fixture(autouse=True)
def clear_db():
    molecules_db.clear()

def test_substructure_search():
    molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"
    expected_result = ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    result = substructure_search(molecules, substructure)
    assert result == expected_result

def test_add_molecule_to_db():
    add_molecule_to_db("mol1", "CCO")
    assert molecules_db["mol1"] == "CCO"
    with pytest.raises(HTTPException):
        add_molecule_to_db("mol1", "CCO")

def test_add_molecule():
    response = client.post("/add_molecule/", json={"identifier": "mol1", "smiles": "CCO"})
    assert response.status_code == 201
    assert response.json() == {"message": "Molecule successfully added in database"}

def test_get_molecule():
    add_molecule_to_db("mol1", "CCO")
    response = client.get("/get_molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}

def test_update_molecule():
    add_molecule_to_db("mol1", "CCO")
    response = client.put("/update_molecule/mol1", json={"smiles": "OCC"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}
    assert molecules_db["mol1"] == "OCC"

def test_delete_molecule():
    add_molecule_to_db("mol1", "CCO")
    response = client.delete("/delete_molecule/mol1")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully"}
    assert "mol1" not in molecules_db

def test_list_molecules():
    add_molecule_to_db("mol1", "CCO")
    add_molecule_to_db("mol2", "OCC")
    response = client.get("/list_molecules/")
    assert response.status_code == 200
    assert response.json() == {"mol1": "CCO", "mol2": "OCC"}

def test_substructure_search_api():
    add_molecule_to_db("mol1", "c1ccccc1")
    add_molecule_to_db("mol2", "CC(=O)Oc1ccccc1C(=O)O")
    response = client.post("/substructure_search/", json={"mol": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"matching_molecules": ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]}

def test_upload_molecules(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "test.csv"
    p.write_text("identifier,smiles\nmol1,CCO\nmol2,OCC\n")
    response = client.post("/upload_molecules/", json={"file_path": str(p)})
    assert response.status_code == 201
    assert response.json() == {"message": "2 molecules from external database uploaded successfully"}
    assert molecules_db["mol1"] == "CCO"
    assert molecules_db["mol2"] == "OCC"
