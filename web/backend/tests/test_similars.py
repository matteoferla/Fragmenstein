"""Test similars endpoints — manual SMILES paste and file upload."""

import io


def _create_session_with_demo(client) -> str:
    r = client.post("/api/sessions", json={"name": "similars-test"})
    sid = r.json()["id"]
    client.post(f"/api/sessions/{sid}/demo/mpro?n_hits=3")
    return sid


def test_manual_smiles(client):
    r = client.post("/api/sessions", json={"name": "manual-test"})
    sid = r.json()["id"]

    smiles_text = "c1ccccc1 benzene\nCC(=O)O acetic_acid\nCCO ethanol\nINVALID_SMILES bad_one"
    r = client.post(
        f"/api/sessions/{sid}/similars/manual",
        json={"smiles_text": smiles_text},
    )
    assert r.status_code == 200
    data = r.json()
    assert data["valid"] == 3
    assert data["invalid"] == 1
    assert "job_id" in data

    # Results should be fetchable
    r = client.get(f"/api/jobs/{data['job_id']}/results")
    assert r.status_code == 200
    assert r.json()["count"] == 3


def test_manual_smiles_empty(client):
    r = client.post("/api/sessions", json={"name": "empty-test"})
    sid = r.json()["id"]

    r = client.post(
        f"/api/sessions/{sid}/similars/manual",
        json={"smiles_text": "   "},
    )
    assert r.status_code == 400


def test_upload_csv(client):
    r = client.post("/api/sessions", json={"name": "csv-test"})
    sid = r.json()["id"]

    csv_content = "smiles,name\nc1ccccc1,benzene\nCCO,ethanol\nCC(=O)O,acetic_acid\n"
    r = client.post(
        f"/api/sessions/{sid}/similars/upload",
        files={"file": ("compounds.csv", io.BytesIO(csv_content.encode()), "text/csv")},
    )
    assert r.status_code == 200
    data = r.json()
    assert data["valid"] == 3
    assert data["invalid"] == 0

    r = client.get(f"/api/jobs/{data['job_id']}/results")
    assert r.json()["count"] == 3


def test_upload_csv_no_smiles_column(client):
    r = client.post("/api/sessions", json={"name": "bad-csv"})
    sid = r.json()["id"]

    csv_content = "molecule,id\nc1ccccc1,1\n"
    r = client.post(
        f"/api/sessions/{sid}/similars/upload",
        files={"file": ("bad.csv", io.BytesIO(csv_content.encode()), "text/csv")},
    )
    assert r.status_code == 400
    assert "SMILES column" in r.json()["detail"]


def test_upload_tsv(client):
    r = client.post("/api/sessions", json={"name": "tsv-test"})
    sid = r.json()["id"]

    tsv_content = "smiles\tname\nc1ccccc1\tbenzene\nCCO\tethanol\n"
    r = client.post(
        f"/api/sessions/{sid}/similars/upload",
        files={"file": ("compounds.tsv", io.BytesIO(tsv_content.encode()), "text/tab-separated-values")},
    )
    assert r.status_code == 200
    assert r.json()["valid"] == 2
