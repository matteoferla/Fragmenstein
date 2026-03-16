"""Test molecule endpoints — depict, molblocks."""


def test_depict_valid_smiles(client):
    r = client.get("/api/depict?smiles=c1ccccc1&width=200&height=150")
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"
    assert len(r.content) > 100  # should be a real PNG


def test_depict_invalid_smiles(client):
    r = client.get("/api/depict?smiles=INVALID_SMILES")
    assert r.status_code == 400
