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
    Batch‑upload all new abstracts as a single JSONL file, embedding metadata
    into the text and using a static chunking strategy so each record is one chunk.
    """
    # 1) Load or init seen‑PMID manifest
    mpath = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(mpath.read_text())) if mpath.exists() else set()

    lines = []
    for row in df.itertuples():
        pmid = row.pmid
        if not pmid or pmid in seen:
            continue

        # Build a single TEXT block including metadata headers
        block = (
            f"DOI:{row.doi or ''}\n"
            f"PMC:{row.pmc or ''}\n"
            f"Title:{row.title}\n"
            "Abstract:\n"
            f"{row.abstract}"
        )
        lines.append({"text": block})
        seen.add(pmid)

    if not lines:
        print(f"→ No new abstracts for {vstore_id}.")
        return

    # 2) Write JSONL
    jsonl = "\n".join(json.dumps(l, ensure_ascii=False) for l in lines)
    buf = io.BytesIO(jsonl.encode("utf-8"))
    buf.name = "abstracts_batch.jsonl"

    # 3) Upload & attach with a static chunking strategy
    #    so each JSONL line becomes exactly one vector chunk
    client.vector_stores.files.upload_and_poll(
        vector_store_id=vstore_id,
        file=buf,
        chunking_strategy={
            "type": "static",
            "static": {
                "max_chunk_size_tokens": 100_000,
                "chunk_overlap_tokens": 0
            }
        }
    )

    # 4) Save manifest
    mpath.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {len(lines)} abstracts in one batch into {vstore_id}")