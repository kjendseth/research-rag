# config.py

import os
import json
from pathlib import Path
from dataclasses import dataclass, field

import openai

@dataclass
class Settings:
    openai_api_key: str = field(default_factory=lambda: os.getenv("OPENAI_API_KEY"))
    ncbi_email: str     = field(default_factory=lambda: os.getenv("NCBI_EMAIL"))
    wos_api_key: str | None = os.getenv("WOS_API_KEY")
    ezproxy_prefix: str | None = os.getenv("EZPROXY_PREFIX")

    # feature toggles
    refine_search: bool = True
    use_headless_browser: bool = False
    embedding_model: str = "text-embedding-3-large"
    llm_model: str = "gpt-4.1"
    chunk_size: int = 800
    chunk_overlap: int = 100

# instantiate settings once
settings = Settings()

def get_store_ids(project: str) -> dict[str,str]:
    """
    Ensure there's a store_ids.json under the project folder containing
    valid OpenAI vector store IDs (beginning with 'vs_'). If missing, create
    two new stores via the API and save them.
    """
    workdir = Path(project)
    workdir.mkdir(exist_ok=True)
    sid_file = workdir / "store_ids.json"

    client = openai.OpenAI(api_key=settings.openai_api_key)

    if sid_file.exists():
        return json.loads(sid_file.read_text())

    # Create the vector stores (name only; description not supported)
    abstracts_vs = client.vector_stores.create(
        name=f"{project}-abstracts"
    ).id
    pdfs_vs = client.vector_stores.create(
        name=f"{project}-fulltext"
    ).id

    ids = {"abstracts": abstracts_vs, "pdfs": pdfs_vs}
    sid_file.write_text(json.dumps(ids, indent=2))
    print(f"✅ Created vector stores: {ids}")
    return ids