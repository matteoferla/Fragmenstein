"""Test upload endpoints — template and hits."""

import io


def _create_session(client) -> str:
    r = client.post("/api/sessions", json={"name": "upload-test"})
    return r.json()["id"]


def test_upload_template(client):
    sid = _create_session(client)
    pdb_content = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n"

    r = client.post(
        f"/api/sessions/{sid}/template",
        files={"file": ("test.pdb", io.BytesIO(pdb_content), "chemical/x-pdb")},
    )
    assert r.status_code == 200
    assert r.json()["filename"] == "test.pdb"

    # Verify template is stored
    r = client.get(f"/api/sessions/{sid}")
    assert r.json()["template_filename"] == "test.pdb"


def test_upload_template_wrong_format(client):
    sid = _create_session(client)

    r = client.post(
        f"/api/sessions/{sid}/template",
        files={"file": ("test.txt", io.BytesIO(b"not a pdb"), "text/plain")},
    )
    assert r.status_code == 400


def test_get_template_pdb(client):
    sid = _create_session(client)
    pdb_text = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n"

    client.post(
        f"/api/sessions/{sid}/template",
        files={"file": ("t.pdb", io.BytesIO(pdb_text.encode()), "chemical/x-pdb")},
    )

    r = client.get(f"/api/sessions/{sid}/template/pdb")
    assert r.status_code == 200
    assert "ATOM" in r.json()["pdb"]


def test_upload_hits_mol(client):
    sid = _create_session(client)

    # Minimal MOL block for methane
    mol_block = """methane
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
    r = client.post(
        f"/api/sessions/{sid}/hits",
        files={"files": ("methane.mol", io.BytesIO(mol_block.encode()), "chemical/x-mdl-molfile")},
    )
    assert r.status_code == 200
    assert "1 molecule" in r.json()["message"]

    # Verify hits list
    r = client.get(f"/api/sessions/{sid}/hits")
    assert r.status_code == 200
    assert r.json()["count"] == 1
    assert r.json()["hits"][0]["name"] == "methane"


def test_list_hits_empty(client):
    sid = _create_session(client)
    r = client.get(f"/api/sessions/{sid}/hits")
    assert r.status_code == 200
    assert r.json()["count"] == 0


def test_clean_template(client):
    sid = _create_session(client)
    pdb_text = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    2  O   HOH A 100       1.000   1.000   1.000  1.00  0.00           O\n"
        "END\n"
    )
    client.post(
        f"/api/sessions/{sid}/template",
        files={"file": ("t.pdb", io.BytesIO(pdb_text.encode()), "chemical/x-pdb")},
    )

    r = client.post(f"/api/sessions/{sid}/template/clean?remove_residues=HOH")
    assert r.status_code == 200
    assert r.json()["removed_count"] == 1

    # Verify HOH is gone
    r = client.get(f"/api/sessions/{sid}/template/pdb")
    assert "HOH" not in r.json()["pdb"]
    assert "ALA" in r.json()["pdb"]
