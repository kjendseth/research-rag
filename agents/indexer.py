# agents/indexer.py

import os
import io
import json
import hashlib
import openai
from pathlib import Path
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)
_ABSTRACT_MANIFEST_PREFIX = ".abstract_manifest_"
_PDF_MANIFEST = ".pdf_manifest.json"

def index_abstracts(df, vstore_id: str):
    """
    Incrementally upload new abstracts (with full metadata) including keywords.
    """
    manifest_path = Path(f"{_ABSTRACT_MANIFEST_PREFIX}{vstore_id}.json")
    seen = set(json.loads(manifest_path.read_text())) if manifest_path.exists() else set()

    new_records = []
    for row in df.itertuples():
        doi = getattr(row, "doi", None)
        if not doi or doi in seen:
            continue
        new_records.append({
            "doi": doi,
            "title": row.title,
            "abstract": row.abstract,
            "authors": row.authors,
            "journal": row.journal,
            "year": row.year,
            "keywords": row.keywords
        })
        seen.add(doi)

    if not new_records:
        print(f"→ No new abstracts to ingest for store {vstore_id}.")
        return

    manifest_path.write_text(json.dumps(sorted(seen)))

    payload = json.dumps(new_records).encode("utf-8")
    file_obj = io.BytesIO(payload)
    file_obj.name = "new_abstracts.json"

    resp = client.files.create(
        file=file_obj,
        purpose="user_data"
    )
    file_id = resp.id

    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=file_id
    )
    print(f"＋ Ingested {len(new_records)} abstracts (with keywords) into {vstore_id} (file {file_id})")