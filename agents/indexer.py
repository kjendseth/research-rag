# agents/indexer.py

import os
import io
import json
import hashlib
import openai
from pathlib import Path
from config import settings

# Initialize OpenAI client
client = openai.OpenAI(api_key=settings.openai_api_key)

# Manifest filenames
_ABSTRACT_MANIFEST_PREFIX = ".abstract_manifest_"
_PDF_MANIFEST = ".pdf_manifest.json"

def index_abstracts(df, vstore_id: str):
    """
    Incrementally upload new abstracts (with full metadata) including keywords.
    Keeps a manifest of seen DOIs under .abstract_manifest_<vstore_id>.json.
    """
    # Load or init DOI manifest
    manifest_path = Path(f"{_ABSTRACT_MANIFEST_PREFIX}{vstore_id}.json")
    if manifest_path.exists():
        seen = set(json.loads(manifest_path.read_text()))
    else:
        seen = set()

    # Collect new records
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

    # Update manifest
    manifest_path.write_text(json.dumps(sorted(seen)))

    # Upload new abstracts as JSON
    payload = json.dumps(new_records).encode("utf-8")
    file_obj = io.BytesIO(payload)
    file_obj.name = "new_abstracts.json"

    resp = client.files.create(
        file=file_obj,
        purpose="user_data"
    )
    file_id = resp.id

    # Register file with vector store
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=file_id
    )
    print(f"＋ Ingested {len(new_records)} abstracts into {vstore_id} (file {file_id})")


def index_pdf(pdf_path: str, vstore_id: str):
    """
    Upload a PDF to the vector store, skipping duplicates by SHA-256.
    Stores the PDF’s DOI (derived from the filename) as metadata.
    """
    path = Path(pdf_path)
    # Compute SHA-256 of file contents
    sha = hashlib.sha256(path.read_bytes()).hexdigest()

    # Load or init PDF manifest
    if Path(_PDF_MANIFEST).exists():
        manifest = json.loads(Path(_PDF_MANIFEST).read_text())
    else:
        manifest = {}

    seen = manifest.get(vstore_id, [])
    if sha in seen:
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # Upload PDF file
    with open(pdf_path, "rb") as f:
        resp = client.files.create(file=f, purpose="user_data")
    file_id = resp.id

    # Attach to vector store with DOI metadata
    doi = path.stem.replace("_", "/")
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=file_id,
        attributes={"doi": doi}
    )
    print(f"＋ Indexed PDF: {pdf_path} (file {file_id})")

    # Update and save manifest
    seen.append(sha)
    manifest[vstore_id] = seen
    Path(_PDF_MANIFEST).write_text(json.dumps(manifest))