"""Job lifecycle management and SSE broadcast."""

import asyncio
import json
import logging
from collections import defaultdict
from datetime import datetime
from typing import AsyncIterator

from ..models.job import Job, create_job, get_job, update_job

log = logging.getLogger(__name__)

# In-memory subscribers for SSE per job
_subscribers: dict[str, list[asyncio.Queue]] = defaultdict(list)


def start_job(session_id: str, job_type: str, config: dict | None = None) -> Job:
    """Create a new job record."""
    return create_job(session_id, job_type, config)


def mark_running(job_id: str):
    update_job(job_id, status="running", started_at=datetime.utcnow().isoformat())


def mark_completed(job_id: str, result_path: str | None = None):
    update_job(
        job_id,
        status="completed",
        progress=1.0,
        message="Done",
        result_path=result_path,
        completed_at=datetime.utcnow().isoformat(),
    )
    _broadcast(job_id, {"progress": 1.0, "status": "completed", "message": "Done"})


def mark_failed(job_id: str, error: str):
    update_job(
        job_id,
        status="failed",
        error=error,
        completed_at=datetime.utcnow().isoformat(),
    )
    _broadcast(job_id, {"progress": 0, "status": "failed", "message": error})


def update_progress(job_id: str, progress: float, message: str = ""):
    update_job(job_id, progress=progress, message=message)
    _broadcast(job_id, {"progress": progress, "status": "running", "message": message})


def _broadcast(job_id: str, data: dict):
    """Push an event to all SSE subscribers for this job."""
    for queue in _subscribers.get(job_id, []):
        try:
            queue.put_nowait(data)
        except asyncio.QueueFull:
            pass  # Drop if consumer is slow


async def subscribe(job_id: str) -> AsyncIterator[str]:
    """Async generator yielding SSE events for a job."""
    queue: asyncio.Queue = asyncio.Queue(maxsize=50)
    _subscribers[job_id].append(queue)
    try:
        # Send current state first
        job = get_job(job_id)
        if job:
            yield f"data: {json.dumps({'progress': job.progress, 'status': job.status, 'message': job.message})}\n\n"
            if job.status in ("completed", "failed"):
                return

        while True:
            try:
                data = await asyncio.wait_for(queue.get(), timeout=30)
                yield f"data: {json.dumps(data)}\n\n"
                if data.get("status") in ("completed", "failed"):
                    return
            except asyncio.TimeoutError:
                # Send keepalive
                yield f": keepalive\n\n"
    finally:
        _subscribers[job_id].remove(queue)
        if not _subscribers[job_id]:
            del _subscribers[job_id]
