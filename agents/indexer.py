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
    Batch‑upload new abstracts as a JSON file where each line has both
    'text' and 'metadata' keys. Static chunking keeps each record intact.
    """
    # 1) load or init seen PMIDs
    mpath = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(mpath.read_text())) if mpath.exists() else set()

    records = []
    for row in df.itertuples():
        pmid = row.pmid
        if not pmid or pmid in seen:
            continue

        # metadata dict
        attrs = {
            "pmid":    pmid,
            "doi":     row.doi or "",
            "pmc":     getattr(row, "pmc", "") or "",
            "title":   row.title or "",
            "journal": row.journal or "",
            "year":    str(row.year or ""),
            "authors": "; ".join(row.authors or []),
            "keywords":"; ".join(row.keywords or [])
        }

        # text payload (can be just abstract if you prefer)
        text_block = f"Title: {row.title}\n\nAbstract:\n{row.abstract}"

        records.append({
            "text":     text_block,
            "metadata": attrs
        })
        seen.add(pmid)

    if not records:
        print(f"→ No new abstracts to ingest for store {vstore_id}.")
        return

    # 2) serialize JSONL and wrap as .json
    jsonl = "\n".join(json.dumps(r, ensure_ascii=False) for r in records)
    buf = io.BytesIO(jsonl.encode("utf-8"))
    buf.name = "abstracts_with_metadata.json"

    # 3) batch import with static chunking (each record → one chunk + attributes)
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

    # 4) save manifest
    mpath.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {len(records)} abstracts in one batch into {vstore_id}")