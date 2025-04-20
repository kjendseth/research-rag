# agents/indexer.py

import io
import json
import hashlib
import openai
from pathlib import Path
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

_ABSTRACT_MANIFEST = ".abstract_manifest.json"
_PDF_MANIFEST      = ".pdf_manifest.json"

def index_abstracts(df, vstore_id: str):
    """
    Incrementally upload each abstract as its own vector‑store file,
    serializing list attributes into strings.
    """
    manifest_path = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(manifest_path.read_text())) if manifest_path.exists() else set()

    new_count = 0
    for row in df.itertuples():
        pmid = getattr(row, "pmid", None)
        if not pmid or pmid in seen:
            continue

        # Serialize list fields into strings
        authors_str  = "; ".join(row.authors or [])
        keywords_str = "; ".join(row.keywords or [])

        attrs = {
            "pmid":    pmid,
            "doi":     row.doi   or "",
            "pmc":     getattr(row, "pmc", "") or "",
            "title":   row.title or "",
            "journal": row.journal or "",
            "year":    str(row.year or ""),
            "authors": authors_str,
            "keywords": keywords_str
        }

        # Build the payload
        payload = json.dumps({"text": f"{row.title}\n\n{row.abstract}"}).encode("utf-8")
        buf = io.BytesIO(payload)
        buf.name = f"{pmid}.json"

        # Upload file
        f = client.files.create(file=buf, purpose="user_data")
        # Attach to vector store with serialized attributes
        client.vector_stores.files.create(
            vector_store_id=vstore_id,
            file_id=f.id,
            attributes=attrs
        )

        seen.add(pmid)
        new_count += 1

    # Save updated manifest
    manifest_path.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {new_count} new abstracts into {vstore_id}")

def index_pdf(pdf_path: str, vstore_id: str):
    """
    Upload a PDF to the vector store, skipping duplicates by SHA‑256,
    attaching DOI in metadata.
    """
    sha = hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest()

    pm = Path(_PDF_MANIFEST)
    manifest = json.loads(pm.read_text()) if pm.exists() else {}

    seen = set(manifest.get(vstore_id, []))
    if sha in seen:
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # Upload PDF
    with open(pdf_path, "rb") as fbin:
        f = client.files.create(file=fbin, purpose="user_data")

    doi = Path(pdf_path).stem.replace("_", "/")
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id,
        attributes={"doi": doi}
    )
    print(f"＋ Indexed PDF {pdf_path} into {vstore_id}")

    # Update PDF manifest
    seen.add(sha)
    manifest[vstore_id] = list(seen)
    pm.write_text(json.dumps(manifest))