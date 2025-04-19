# config.py

from dataclasses import dataclass, field
import os
import json
from pathlib import Path
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
    llm_model: str = "gpt-4o-mini"
    chunk_size: int = 800
    chunk_overlap: int = 100

settings = Settings()

# ——————————————————————————————————————————————————————————————————————————
# Project‐scoped vector store creation & persistence

def get_store_ids(project: str) -> dict[str,str]:
    """
    Ensure there's a store_ids.json under the project folder containing
    real OpenAI vector store IDs (beginning with 'vs_'). If missing, create
    two new stores via the API and save them.
    """
    # project workspace
    workdir = Path(project)
    workdir.mkdir(exist_ok=True)
    sid_file = workdir / "store_ids.json"

    client = openai.OpenAI(api_key=settings.openai_api_key)

    if sid_file.exists():
        return json.loads(sid_file.read_text())

    # Create two vector stores: abstracts & pdfs
    abstracts_vs = client.vector_stores.create(
        name=f"{project}-abstracts", 
        description="Abstract embeddings store"
    ).id
    pdfs_vs = client.vector_stores.create(
        name=f"{project}-fulltext", 
        description="Full‑text PDF embeddings store"
    ).id

    ids = {"abstracts": abstracts_vs, "pdfs": pdfs_vs}
    sid_file.write_text(json.dumps(ids, indent=2))
    print(f"✅ Created vector stores: {ids}")
    return ids