"""Application configuration using Pydantic settings."""

import os
from pathlib import Path
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    app_name: str = "Fragmenstein Web UI"
    debug: bool = True

    # Paths
    data_dir: Path = Path(os.environ.get("FRAG_DATA_DIR", "./data"))
    upload_dir: Path = Path(os.environ.get("FRAG_UPLOAD_DIR", "./data/uploads"))
    work_dir: Path = Path(os.environ.get("FRAG_WORK_DIR", "./data/work"))
    db_path: Path = Path(os.environ.get("FRAG_DB_PATH", "./data/fragmenstein.db"))

    # Limits
    max_upload_size_mb: int = 100
    max_hits: int = 50
    default_n_cores: int = -1
    default_timeout: int = 240
    default_combination_size: int = 2

    # CORS
    cors_origins: list[str] = ["http://localhost:3000"]

    # Victor defaults
    default_victor_type: str = "Wictor"
    default_joining_cutoff: float = 5.0

    model_config = {"env_prefix": "FRAG_"}

    def ensure_dirs(self):
        """Create required directories if they don't exist."""
        for d in [self.data_dir, self.upload_dir, self.work_dir]:
            d.mkdir(parents=True, exist_ok=True)


settings = Settings()
