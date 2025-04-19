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

def index_pdf(pdf_path: str, vstore_id: str):
    """
    Read the PDF, chunk it, dedupe by SHA256, embed and store.
    Skips if the same SHA is already in the manifest for this store,
    or if PDF parsing fails.
    """
    # Compute hash
    with open(pdf_path, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()

    manifest = _load_manifest()
    if sha in manifest.get(vstore_id, []):
        print(f"→ Skipping duplicate PDF: {pdf_path}")
        return

    # Attempt to open PDF
    try:
        reader = PdfReader(pdf_path)
    except (pdf_errors.PdfReadError, pdf_errors.PdfStreamError) as e:
        print(f"⚠️  Skipping invalid PDF {pdf_path}: {e}")
        return

    # Extract text safely
    text_pages = []
    for page in reader.pages:
        try:
            txt = page.extract_text() or ""
            text_pages.append(txt)
        except Exception:
            continue
    full_text = "\n".join(text_pages)

    # If there's no text, skip
    if not full_text.strip():
        print(f"⚠️  No extractable text in {pdf_path}, skipping.")
        return

    # Chunk & embed
    chunks = chunk_text(full_text)
    embeddings = embed(chunks)

    # Send to vector store
    client.vector_stores.create_chunk_embeddings(
        vstore_id,
        embeddings=embeddings,
        metadata=[{"pdf": pdf_path, "sha256": sha} for _ in embeddings],
    )

    # Update manifest
    manifest.setdefault(vstore_id, []).append(sha)
    _save_manifest(manifest)
    print(f"＋ Indexed PDF: {pdf_path} ({len(chunks)} chunks)")
