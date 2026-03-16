"""Shared test fixtures."""

import os
import tempfile

import pytest

# Use a temp directory for test data
_tmpdir = tempfile.mkdtemp(prefix="frag_test_")
os.environ["FRAG_DATA_DIR"] = _tmpdir
os.environ["FRAG_UPLOAD_DIR"] = os.path.join(_tmpdir, "uploads")
os.environ["FRAG_WORK_DIR"] = os.path.join(_tmpdir, "work")
os.environ["FRAG_DB_PATH"] = os.path.join(_tmpdir, "test.db")

from app.main import app  # noqa: E402
from app.database import init_db  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def setup_db():
    """Initialize the test database once."""
    from app.config import settings
    settings.ensure_dirs()
    init_db()


@pytest.fixture
def client():
    """FastAPI TestClient for synchronous tests."""
    from starlette.testclient import TestClient
    with TestClient(app) as c:
        yield c
