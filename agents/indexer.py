# indexer.py  
import json, hashlib, os, pypdf, openai, tiktoken
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)
_tkm = tiktoken.get_encoding("cl100k_base")

_MANIFEST = ".pdf_manifest.json"   # local cache of sha256 → store_ids

def _load_manifest():
    if os.path.exists(_MANIFEST):
        return json.load(open(_MANIFEST))
    return {}

def _save_manifest(m):
    json.dump(m, open(_MANIFEST,"w"))

def index_pdf(pdf_path: str, vstore_id: str):
    sha = hashlib.sha256(open(pdf_path,"rb").read()).hexdigest()
    manifest = _load_manifest()
    if sha in manifest.get(vstore_id, []):
        print(f"✓ duplicate PDF, already in store: {pdf_path}")
        return
    text = "\n".join(p.extract_text() for p in pypdf.PdfReader(pdf_path).pages)
    chunks = chunk_text(text)
    embeds = embed(chunks)
    client.vector_stores.create_chunk_embeddings(
        vstore_id, embeddings=embeds,
        metadata={"pdf": pdf_path, "sha256": sha})
    manifest.setdefault(vstore_id, []).append(sha)
    _save_manifest(manifest)
    print(f"＋ PDF embedded: {pdf_path}")
