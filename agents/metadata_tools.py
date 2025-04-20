# agents/metadata_tools.py

import openai
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

def semantic_search(vector_store_id: str, query: str, limit: int = 100) -> list[dict]:
    """
    Return up to `limit` papers whose abstract/title matches `query` semantically.
    Each dict has keys: doi, pmid, title, snippet.
    """
    # 1) embed the query
    emb = client.embeddings.create(
        model=settings.embedding_model,
        input=[query]
    ).data[0].embedding

    # 2) vector search
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query_embedding=emb,
        limit=limit
    )

    results = []
    for hit in resp.data:
        meta = hit.metadata
        text = hit.document or ""
        snippet = text[:200].replace("\n"," ")
        results.append({
            "doi": meta.get("doi"),
            "pmid": meta.get("pmid"),
            "title": meta.get("title"),
            "snippet": snippet
        })
    return results