"""Test job endpoints — status, cancel."""

from app.services import job_manager


def test_job_status(client):
    r = client.post("/api/sessions", json={"name": "job-test"})
    sid = r.json()["id"]

    job = job_manager.start_job(sid, "combine")
    r = client.get(f"/api/jobs/{job.id}")
    assert r.status_code == 200
    assert r.json()["status"] == "pending"
    assert r.json()["type"] == "combine"


def test_job_not_found(client):
    r = client.get("/api/jobs/nonexistent")
    assert r.status_code == 404


def test_cancel_job(client):
    r = client.post("/api/sessions", json={"name": "cancel-test"})
    sid = r.json()["id"]

    job = job_manager.start_job(sid, "combine")
    job_manager.mark_running(job.id)

    r = client.post(f"/api/jobs/{job.id}/cancel")
    assert r.status_code == 200
    assert r.json()["status"] == "cancelled"

    # Verify status persisted
    r = client.get(f"/api/jobs/{job.id}")
    assert r.json()["status"] == "cancelled"


def test_cancel_completed_job(client):
    r = client.post("/api/sessions", json={"name": "done-test"})
    sid = r.json()["id"]

    job = job_manager.start_job(sid, "combine")
    job_manager.mark_completed(job.id)

    r = client.post(f"/api/jobs/{job.id}/cancel")
    assert r.status_code == 200
    assert r.json()["status"] == "completed"  # unchanged
