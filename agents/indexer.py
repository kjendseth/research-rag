# agents/indexer.py

import io
import json
import openai
from pathlib import Path
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

_ABSTRACT_MANIFEST = ".abstract_manifest.json"
_PDF_MANIFEST      = ".pdf_manifest.json"
_MAX_ATTR_LEN      = 512

def _truncate(s: str) -> str:
    return s if len(s) <= _MAX_ATTR_LEN else s[:_MAX_ATTR_LEN]

def index_abstracts(df, vstore_id: str):
    """
    Batch‑upload new abstracts as a single JSONL file with embedded metadata.
    """
    # 1. load or init seen manifest
    mpath = Path(_ABSTRACT_MANIFEST)
    seen = set(json.loads(mpath.read_text())) if mpath.exists() else set()

    # 2. build list of new-record dicts
    new_meta = []
    for row in df.itertuples():
        pmid = getattr(row, "pmid", None)
        if not pmid or pmid in seen:
            continue

        # prepare metadata, serializing lists to semicolon strings
        attrs = {
            "pmid":    pmid,
            "doi":     row.doi or "",
            "pmc":     getattr(row, "pmc", "") or "",
            "title":   _truncate(row.title or ""),
            "journal": _truncate(row.journal or ""),
            "year":    str(row.year or ""),
            "authors": _truncate("; ".join(row.authors or [])),
            "keywords":_truncate("; ".join(row.keywords or []))
        }

        new_meta.append({
            "text": f"{row.title}\n\n{row.abstract}",
            "metadata": attrs
        })
        seen.add(pmid)

    if not new_meta:
        print(f"→ No new abstracts to ingest for store {vstore_id}.")
        return

    # 3. serialize to JSONL
    jsonl = "\n".join(json.dumps(rec, ensure_ascii=False) for rec in new_meta)
    buf = io.BytesIO(jsonl.encode("utf-8"))
    buf.name = "new_abstracts.jsonl"

    # 4. upload file once
    f = client.files.create(file=buf, purpose="user_data")
    # 5. import entire batch into vector store
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id
    )

    # 6. write updated manifest
    mpath.write_text(json.dumps(sorted(seen)))
    print(f"＋ Indexed {len(new_meta)} abstracts in one batch into {vstore_id}")

def index_pdf(pdf_path: str, vstore_id: str):
    """
    (Unchanged) Upload a single PDF if not already seen.
    """
    import hashlib
    sha = hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest()

    pm = Path(_PDF_MANIFEST)
    manifest = json.loads(pm.read_text()) if pm.exists() else {}
    seen = set(manifest.get(vstore_id, []))

    if sha in seen:
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    with open(pdf_path, "rb") as fbin:
        f = client.files.create(file=fbin, purpose="user_data")

    doi = Path(pdf_path).stem.replace("_", "/")
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=f.id,
        attributes={"doi": doi}
    )
    print(f"＋ Indexed PDF {pdf_path} into {vstore_id}")

    seen.add(sha)
    manifest[vstore_id] = list(seen)
    pm.write_text(json.dumps(manifest))