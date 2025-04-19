# agents/indexer.py

import os
import json
import hashlib
import openai
import io
#import pypdf
from pypdf import PdfReader, errors as pdf_errors
import tiktoken
from tqdm import tqdm
from config import settings

# Initialize OpenAI client
client = openai.OpenAI(api_key=settings.openai_api_key)

# Tokenizer for chunking
_tkm = tiktoken.get_encoding("cl100k_base")

# Local manifest to dedupe PDF embeddings by hash
_MANIFEST = ".pdf_manifest.json"

def _load_manifest():
    if os.path.exists(_MANIFEST):
        return json.load(open(_MANIFEST))
    return {}

def _save_manifest(m):
    json.dump(m, open(_MANIFEST, "w"))

def chunk_text(text: str) -> list[str]:
    """Split large text into overlapping chunks."""
    toks = _tkm.encode(text)
    size, ov = settings.chunk_size, settings.chunk_overlap
    return [
        _tkm.decode(toks[i : i + size])
        for i in range(0, len(toks), size - ov)
    ]

def embed(chunks: list[str]) -> list[list[float]]:
    """Call OpenAI to get embeddings for each chunk."""
    resp = client.embeddings.create(
        input=chunks, model=settings.embedding_model
    )
    return [d.embedding for d in resp.data]

def index_abstracts(df, vstore_id: str):
    """
    Upload all (doi, abstract) records as one JSON file (as user_data),
    then register it with the abstract vector store.
    """
    data = [
        {"doi": row.doi, "abstract": row.abstract}
        for row in df.itertuples()
        if getattr(row, "doi", None)
    ]
    raw = json.dumps(data).encode("utf-8")
    file_obj = io.BytesIO(raw)
    file_obj.name = "abstracts.json"

    # Upload as user_data
    resp = client.files.create(
        file=file_obj,
        purpose="user_data"      
    )
    file_id = resp.id

    # Now link the file into your vector store
    client.vector_stores.files.create(vstore_id, file_id=file_id)
    print(f"✓ Uploaded {len(data)} abstracts to store {vstore_id} (file_id={file_id})")
        
client = openai.OpenAI(api_key=settings.openai_api_key)
_tkm = tiktoken.get_encoding("cl100k_base")
_MANIFEST = ".pdf_manifest.json"

def _load_manifest():
    if os.path.exists(_MANIFEST):
        return json.load(open(_MANIFEST))
    return {}

def _save_manifest(m):
    json.dump(m, open(_MANIFEST, "w"))

def chunk_text(text: str) -> list[str]:
    toks = _tkm.encode(text)
    size, ov = settings.chunk_size, settings.chunk_overlap
    return [
        _tkm.decode(toks[i : i + size])
        for i in range(0, len(toks), size - ov)
    ]

def embed(chunks: list[str]) -> list[list[float]]:
    resp = client.embeddings.create(
        input=chunks, model=settings.embedding_model
    )
    return [d.embedding for d in resp.data]

client = openai.OpenAI(api_key=settings.openai_api_key)
_MANIFEST = ".pdf_manifest.json"

def _load_manifest():
    if os.path.exists(_MANIFEST):
        import json
        return json.load(open(_MANIFEST))
    return {}

def _save_manifest(m):
    import json
    json.dump(m, open(_MANIFEST, "w"))

def index_pdf(pdf_path: str, vstore_id: str):
    """
    Upload the PDF file to OpenAI (purpose=user_data) and register it
    with the given vector store. The API will do chunking+embedding internally.
    """
    # Skip duplicates by SHA‑256
    sha = hashlib.sha256(open(pdf_path, "rb").read()).hexdigest()
    manifest = _load_manifest()
    if sha in manifest.get(vstore_id, []):
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # Upload PDF to file storage
    with open(pdf_path, "rb") as f:
        resp = client.files.create(file=f, purpose="user_data")
    file_id = resp.id

    # Attach file to vector store
    client.vector_stores.files.create(
        vector_store_id=vstore_id,
        file_id=file_id,
    )
    print(f"＋ Indexed PDF: {pdf_path} as file {file_id}")

    # Update manifest
    manifest.setdefault(vstore_id, []).append(sha)
    _save_manifest(manifest)
