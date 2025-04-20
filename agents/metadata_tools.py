# agents/metadata_tools.py

import openai
from config import settings

client = openai.OpenAI(api_key=settings.openai_api_key)

def semantic_search(vector_store_id: str, query: str, limit: int = 100) -> list[dict]:
    """
    Perform a semantic search (via `query` string) and then extract
    DOI, Title, and one-line snippet from the chunk text itself.
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query=query
    )
    hits = resp.data[:limit]
    results = []

    for hit in hits:
        # Get the raw chunk text (fallback from document to text)
        raw = getattr(hit, "document", None) or getattr(hit, "text", "") or ""
        lines = raw.splitlines()

        # Parse DOI from the first line "DOI:..."
        doi = ""
        if lines and lines[0].startswith("DOI:"):
            doi = lines[0][4:].strip()

        # Parse Title from the line that starts with "Title:"
        title = ""
        for line in lines:
            if line.startswith("Title:"):
                title = line[6:].strip()
                break

        # Parse Abstract snippet: first sentence after "Abstract:"
        snippet = ""
        if "Abstract:" in raw:
            abstract_part = raw.split("Abstract:", 1)[1].strip()
            # take up to first period
            snippet = abstract_part.split(".", 1)[0].strip() + "."

        results.append({
            "doi": doi,
            "title": title,
            "snippet": snippet
        })

    return results

def search_by_author(vector_store_id: str, author_substr: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries whose inline chunk text contains the author substring.
    (Assumes you embedded 'Title:' and 'Authors:' lines in your text blocks.)
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query=author_substr
    )
    hits = resp.data[:limit]
    return [{"doi": d["doi"], "title": d["title"], "snippet": d["snippet"]}
            for d in semantic_search(vector_store_id, author_substr, limit)]

def search_by_pmc(vector_store_id: str, limit: int = 1000) -> list[dict]:
    """
    Return up to `limit` entries whose chunk text contains "PMC".
    """
    resp = client.vector_stores.search(
        vector_store_id=vector_store_id,
        query="PMC"
    )
    hits = resp.data[:limit]
    return [{"doi": d["doi"], "title": d["title"], "snippet": d["snippet"]}
            for d in semantic_search(vector_store_id, "PMC", limit)]