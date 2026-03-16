"""FastAPI application entry point."""

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .database import init_db
from .routers import combine, demo, jobs, molecules, monster, place, sessions, similars, single_victor, upload

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup/shutdown events."""
    settings.ensure_dirs()
    init_db()
    log.info("Fragmenstein Web API started")
    yield
    log.info("Fragmenstein Web API shutting down")


app = FastAPI(
    title=settings.app_name,
    version="0.1.0",
    lifespan=lifespan,
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Routers
app.include_router(sessions.router)
app.include_router(upload.router)
app.include_router(demo.router)
app.include_router(combine.router)
app.include_router(similars.router)
app.include_router(place.router)
app.include_router(monster.router)
app.include_router(single_victor.router)
app.include_router(jobs.router)
app.include_router(molecules.router)


@app.get("/api/health")
def health():
    return {"status": "ok", "app": settings.app_name}


@app.get("/api/system-info")
def system_info():
    import os
    import sys
    import platform
    from pathlib import Path

    # Read version from shared VERSION file
    version_file = Path(__file__).resolve().parent.parent.parent / "VERSION"
    version = version_file.read_text().strip() if version_file.exists() else "0.0.0"

    info: dict = {
        "version": version,
        "python": sys.version.split()[0],
        "platform": platform.machine(),
        "cores": os.cpu_count(),
    }

    # PyRosetta
    try:
        import pyrosetta
        info["pyrosetta"] = True
        info["pyrosetta_version"] = pyrosetta.rosetta.utility.Version.version()
    except Exception:
        info["pyrosetta"] = False

    # RDKit
    try:
        from rdkit import rdBase
        info["rdkit"] = rdBase.rdkitVersion
    except Exception:
        info["rdkit"] = None

    # GPU (CUDA via torch)
    try:
        import torch
        if torch.cuda.is_available():
            info["gpu"] = torch.cuda.get_device_name(0)
    except Exception:
        pass

    return info
