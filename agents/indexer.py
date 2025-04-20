# agents/indexer.py

import io
import json
import openai
from pathlib import Path
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

_ABSTRACT_MANIFEST = ".abstract_manifest.json"

def index_abstracts(df, vstore_id: str):
    """
    Batch‑upload all new abstracts as a single JSON file with inline metadata
    and a static chunking strategy (max 4096 tokens per chunk).
    """
    # 1) Load or initialize seen‑PMID manifest
    mpath = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(mpath.read_text())) if mpath.exists() else set()

    # 2) Build list of new records
    records = []
    for row in df.itertuples():
        pmid = row.pmid
        if not pmid or pmid in seen:
            continue

        block = (
            f"DOI:{row.doi or ''}\n"
            f"PMC:{row.pmc or ''}\n"
            f"Title:{row.title}\n"
            "Abstract:\n"
            f"{row.abstract}"
        )
        records.append({"text": block})
        seen.add(pmid)

    if not records:
        print(f"→ No new abstracts to ingest for store {vstore_id}.")
        return

    # 3) Serialize to JSONL and wrap as .json
    jsonl = "\n".join(json.dumps(r, ensure_ascii=False) for r in records)
    buf = io.BytesIO(jsonl.encode("utf-8"))
    buf.name = "abstracts_batch.json"   # must be .json

    # 4) Upload + import with static chunking (max 4096 tokens per chunk)
    client.vector_stores.files.upload_and_poll(
        vector_store_id=vstore_id,
        file=buf,
        chunking_strategy={
            "type": "static",
            "static": {
                "max_chunk_size_tokens": 4096,
                "chunk_overlap_tokens": 0
            }
        }
    )

    # 5) Save updated manifest
    mpath.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {len(records)} abstracts in one batch into {vstore_id}")