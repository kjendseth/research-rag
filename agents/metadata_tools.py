# agents/metadata_tools.py

import openai
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

def semantic_search(vector_store_id: str, query: str, limit: int = 100) -> list[dict]:
    """
    Return up to `limit` papers whose abstract/title matches `query` semantically.
    Each dict has keys: doi, pmid, title, snippet.
    """
    # 1) call the SDK with the query string
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query=query
    )

    hits = resp.data[:limit]  # enforce limit client‑side
    results = []
    for hit in hits:
        meta = hit.metadata or {}
        text = hit.document or ""
        snippet = text[:200].replace("\n", " ")
        results.append({
            "doi":     meta.get("doi"),
            "pmid":    meta.get("pmid"),
            "title":   meta.get("title"),
            "snippet": snippet
        })
    return results

def search_by_author(vector_store_id: str, author_substr: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries whose metadata.authors list contains `author_substr`.
    """
    filt = {"metadata.authors": {"$contains": author_substr.lower()}}
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query="",       # required by signature
        filter=filt
    )

    hits = resp.data[:limit]
    return [
        {
            "doi":   hit.metadata.get("doi"),
            "pmid":  hit.metadata.get("pmid"),
            "title": hit.metadata.get("title")
        }
        for hit in hits
    ]

def search_by_pmc(vector_store_id: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries that have a non‑empty PMC metadata field.
    """
    filt = {"metadata.pmc": {"$ne": ""}}
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query="",      # required by signature
        filter=filt
    )

    hits = resp.data[:limit]
    return [
        {
            "doi":   hit.metadata.get("doi"),
            "pmid":  hit.metadata.get("pmid"),
            "title": hit.metadata.get("title")
        }
        for hit in hits
    ]