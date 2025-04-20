# agents/metadata_tools.py

import openai
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

def semantic_search(vector_store_id: str, query: str, limit: int = 100) -> list[dict]:
    """
    Return up to `limit` papers whose abstract/title matches `query` semantically.
    Each dict has keys: doi, pmid, title, snippet.
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query=query
    )
    hits = resp.data[:limit]
    results = []
    for hit in hits:
        meta = getattr(hit, "attributes", {}) or {}
        # snippet may be under `hit.document` or `hit.text`
        raw = getattr(hit, "document", None) or getattr(hit, "text", "") or ""
        snippet = raw[:200].replace("\n", " ")
        results.append({
            "doi":     meta.get("doi"),
            "pmid":    meta.get("pmid"),
            "title":   meta.get("title"),
            "snippet": snippet
        })
    return results

def search_by_author(vector_store_id: str, author_substr: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries whose attributes.authors list contains `author_substr`.
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query="",
        filter={"attributes.authors": {"$contains": author_substr.lower()}}
    )
    hits = resp.data[:limit]
    return [
        {
            "doi":   getattr(hit, "attributes", {}).get("doi"),
            "pmid":  getattr(hit, "attributes", {}).get("pmid"),
            "title": getattr(hit, "attributes", {}).get("title")
        }
        for hit in hits
    ]

def search_by_pmc(vector_store_id: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries that have a non‑empty PMC in attributes.
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query="",
        filter={"attributes.pmc": {"$ne": ""}}
    )
    hits = resp.data[:limit]
    return [
        {
            "doi":   getattr(hit, "attributes", {}).get("doi"),
            "pmid":  getattr(hit, "attributes", {}).get("pmid"),
            "title": getattr(hit, "attributes", {}).get("title")
        }
        for hit in hits
    ]