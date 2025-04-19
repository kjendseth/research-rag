# config.py
from dataclasses import dataclass, field
import os, uuid, pathlib, json

@dataclass
class Settings:
    openai_api_key: str = field(default_factory=lambda: os.getenv("OPENAI_API_KEY"))
    ncbi_email: str = field(default_factory=lambda: os.getenv("NCBI_EMAIL"))
    wos_api_key: str | None = os.getenv("WOS_API_KEY")
    ezproxy_prefix: str | None = os.getenv("EZPROXY_PREFIX")

    # ---- feature toggles ----
    refine_search: bool = True          # <‑‑ turn query‑improver on/off
    use_headless_browser: bool = False  # set True if institutional PDFs needed
    embedding_model: str = "text-embedding-3-large"
    llm_model: str = "gpt-4o-mini"
    chunk_size: int = 800    # tokens
    chunk_overlap: int = 100

settings = Settings()


def get_store_ids(project: str):
    pth = pathlib.Path(f"{project}/store_ids.json")
    if pth.exists():
        return json.load(open(pth))
    pth.parent.mkdir(exist_ok=True)
    ids = {
        "abstracts": uuid.uuid4().hex,
        "pdfs": uuid.uuid4().hex
    }
    json.dump(ids, open(pth,"w"))
    return ids
