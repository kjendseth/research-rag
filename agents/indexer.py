# agents/indexer.py

import os
import io
import json
import hashlib
import openai
from pathlib import Path
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

# keep a manifest of seen PMIDs so we ingest incrementally
_ABSTRACT_MANIFEST = ".abstract_manifest.json"
_PDF_MANIFEST      = ".pdf_manifest.json"

def index_abstracts(df, vstore_id: str):
    """
    Incrementally upload each abstract as its own vector-store file,
    attaching metadata attributes so searches return DOI/title/etc.
    """
    # load or init manifest of PMIDs
    mpath = Path(_ABSTRACT_MANIFEST)
    if mpath.exists():
        seen = set(json.loads(mpath.read_text()))
    else:
        seen = set()

    new_count = 0
    for row in df.itertuples():
        pmid = getattr(row, "pmid", None)
        if not pmid or pmid in seen:
            continue

        # build metadata attributes
        attrs = {
            "pmid":    pmid,
            "doi":     row.doi or "",
            "pmc":     getattr(row, "pmc", "") or "",
            "title":   row.title or "",
            "journal": row.journal or "",
            "year":    str(row.year or ""),
            "authors": row.authors or [],
            "keywords": row.keywords or []
        }

        # prepare the text payload
        payload = json.dumps({
            "text": f"{row.title}\n\n{row.abstract}"
        }).encode("utf-8")
        buf = io.BytesIO(payload)
        buf.name = f"{pmid}.json"

        # upload file
        f = client.files.create(file=buf, purpose="user_data")

        # attach to vector store with metadata
        client.vector_stores.files.create(
            vector_store_id=vstore_id,
            file_id=f.id,
            attributes=attrs
        )

        seen.add(pmid)
        new_count += 1

    # update manifest
    mpath.write_text(json.dumps(sorted(seen)))

    print(f"＋ Indexed {new_count} new abstracts into {vstore_id}")

def index_pdf(pdf_path: str, vstore_id: str):
    """
    Upload a PDF to the vector store, skipping duplicates by SHA‑256.
    Stores each PDF’s DOI (from filename) as metadata.
    """
    sha = hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest()

    # load or init PDF manifest
    pm = Path(_PDF_MANIFEST)
    manifest = json.loads(pm.read_text()) if pm.exists() else {}

    seen = set(manifest.get(vstore_id, []))
    if sha in seen:
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # upload PDF
    with open(pdf_path, "rb") as fbin:
        f = client.files.create(file=fbin, purpose="user_data")

    doi = Path(pdf_path).stem.replace("_", "/")
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id,
        attributes={"doi": doi}
    )
    print(f"＋ Indexed PDF {pdf_path} into {vstore_id}")

    # update manifest
    seen.add(sha)
    manifest[vstore_id] = list(seen)
    pm.write_text(json.dumps(manifest))