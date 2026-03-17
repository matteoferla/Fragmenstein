"""Test session CRUD endpoints."""


def test_create_session(client):
    r = client.post("/api/sessions", json={"name": "test-session"})
    assert r.status_code == 201
    data = r.json()
    assert data["name"] == "test-session"
    assert data["status"] == "created"
    assert data["hit_count"] == 0
    assert "id" in data


def test_list_sessions(client):
    # Create two sessions
    client.post("/api/sessions", json={"name": "s1"})
    client.post("/api/sessions", json={"name": "s2"})

    r = client.get("/api/sessions")
    assert r.status_code == 200
    sessions = r.json()
    assert len(sessions) >= 2


def test_get_session(client):
    r = client.post("/api/sessions", json={"name": "get-me"})
    sid = r.json()["id"]

    r = client.get(f"/api/sessions/{sid}")
    assert r.status_code == 200
    assert r.json()["name"] == "get-me"


def test_get_session_not_found(client):
    r = client.get("/api/sessions/nonexistent-id")
    assert r.status_code == 404


def test_delete_session(client):
    r = client.post("/api/sessions", json={"name": "delete-me"})
    sid = r.json()["id"]

    r = client.delete(f"/api/sessions/{sid}")
    assert r.status_code == 204

    r = client.get(f"/api/sessions/{sid}")
    assert r.status_code == 404


def test_list_session_jobs(client):
    r = client.post("/api/sessions", json={"name": "jobs-test"})
    sid = r.json()["id"]

    r = client.get(f"/api/sessions/{sid}/jobs")
    assert r.status_code == 200
    assert r.json() == []
