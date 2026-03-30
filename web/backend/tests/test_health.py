"""Test health and system-info endpoints."""


def test_health(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    data = r.json()
    assert data["status"] == "ok"


def test_system_info(client):
    r = client.get("/api/system-info")
    assert r.status_code == 200
    data = r.json()
    assert "version" in data
    assert "python" in data
    assert "cores" in data
    assert isinstance(data["cores"], int)
    assert data["cores"] > 0
