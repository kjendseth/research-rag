# agents/indexer.py

import io
import json
import hashlib
import openai
from pathlib import Path
from config import settings

# Initialize OpenAI client
client = openai.OpenAI(api_key=settings.openai_api_key)

# Manifest filenames
_ABSTRACT_MANIFEST = ".abstract_manifest.json"
_PDF_MANIFEST      = ".pdf_manifest.json"
_MAX_ATTR_LEN      = 512

def _truncate(s: str) -> str:
    """Ensure attribute string ≤ 512 chars."""
    return s if len(s) <= _MAX_ATTR_LEN else s[:_MAX_ATTR_LEN]

def index_abstracts(df, vstore_id: str):
    """
    Batch‑upload new abstracts as a single JSON file with embedded metadata.
    """
    # 1) Load or initialize seen‐PMID manifest
    mpath = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(mpath.read_text())) if mpath.exists() else set()

    # 2) Build list of new records
    new_meta = []
    for row in df.itertuples():
        pmid = getattr(row, "pmid", None)
        if not pmid or pmid in seen:
            continue

        # Serialize and truncate list fields
        authors_str  = _truncate("; ".join(row.authors or []))
        keywords_str = _truncate("; ".join(row.keywords or []))

        # Prepare metadata attributes
        attrs = {
            "pmid":     pmid,
            "doi":      row.doi   or "",
            "pmc":      getattr(row, "pmc", "") or "",
            "title":    _truncate(row.title or ""),
            "journal":  _truncate(row.journal or ""),
            "year":     str(row.year or ""),
            "authors":  authors_str,
            "keywords": keywords_str
        }

        # Append the record
        new_meta.append({
            "text":     f"{row.title}\n\n{row.abstract}",
            "metadata": attrs
        })
        seen.add(pmid)

    if not new_meta:
        print(f"→ No new abstracts to ingest for store {vstore_id}.")
        return

    # 3) Serialize to JSONL and wrap as a .json file
    jsonl = "\n".join(json.dumps(rec, ensure_ascii=False) for rec in new_meta)
    buf = io.BytesIO(jsonl.encode("utf-8"))
    buf.name = "new_abstracts.json"  # must use .json extension

    # 4) Upload the batch file
    f = client.files.create(file=buf, purpose="user_data")

    # 5) Import entire batch into vector store
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id
    )

    # 6) Update and save manifest
    mpath.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {len(new_meta)} abstracts in one batch into {vstore_id}")

def index_pdf(pdf_path: str, vstore_id: str):
    """
    Upload a PDF to the vector store, skipping duplicates by SHA‑256.
    Stores each PDF’s DOI (from filename) as metadata.
    """
    path = Path(pdf_path)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()

    # Load or initialize PDF manifest
    pm = Path(_PDF_MANIFEST)
    manifest = json.loads(pm.read_text()) if pm.exists() else {}
    seen = set(manifest.get(vstore_id, []))

    if sha in seen:
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # Upload PDF file
    with open(pdf_path, "rb") as fbin:
        f = client.files.create(file=fbin, purpose="user_data")

    # Attach to vector store with DOI metadata
    doi = path.stem.replace("_", "/")
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id,
        attributes={"doi": doi}
    )
    print(f"＋ Indexed PDF {pdf_path} into {vstore_id}")

    # Update and save manifest
    seen.add(sha)
    manifest[vstore_id] = list(seen)
    pm.write_text(json.dumps(manifest))