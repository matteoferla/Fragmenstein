"""Test demo dataset loading."""


def test_load_mpro_demo(client):
    r = client.post("/api/sessions", json={"name": "demo-test"})
    sid = r.json()["id"]

    r = client.post(f"/api/sessions/{sid}/demo/mpro?n_hits=3")
    assert r.status_code == 200
    data = r.json()
    assert data["hit_count"] == 3
    assert "mpro" in data["template_filename"].lower()

    # Verify hits loaded
    r = client.get(f"/api/sessions/{sid}/hits")
    assert r.json()["count"] == 3

    # Verify template loaded
    r = client.get(f"/api/sessions/{sid}/template/pdb")
    assert r.status_code == 200
    assert "ATOM" in r.json()["pdb"]


def test_load_mac1_demo(client):
    r = client.post("/api/sessions", json={"name": "mac1-test"})
    sid = r.json()["id"]

    r = client.post(f"/api/sessions/{sid}/demo/mac1?n_hits=2")
    assert r.status_code == 200
    assert r.json()["hit_count"] == 2


def test_load_unknown_demo(client):
    r = client.post("/api/sessions", json={"name": "bad-demo"})
    sid = r.json()["id"]

    r = client.post(f"/api/sessions/{sid}/demo/unknown")
    assert r.status_code == 400
