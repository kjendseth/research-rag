import requests, pandas as pd, time, hashlib, json, openai, logging, os
from Bio import Entrez          # pip install biopython
from config import settings

_XREF_HEADERS = {
    "User-Agent": f"research-rag/0.1 (mailto:{os.getenv('NCBI_EMAIL','user@example.com')})",
    "Accept": "application/json"
}


def pubmed_search(query, retmax=200):
    Entrez.email = settings.ncbi_email
    ids = Entrez.read(Entrez.esearch(db="pubmed", term=query,
                                     retmax=retmax))["IdList"]
    fetch = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    articles = Entrez.read(fetch)["PubmedArticle"]
    rows = []
    for art in articles:
        art = art["MedlineCitation"]["Article"]
        doi = next((i for i in art["ELocationID"] if i.attributes["EIdType"]=="doi"), None)
        rows.append({
            "title": art["ArticleTitle"],
            "abstract": art.get("Abstract", {}).get("AbstractText", [""])[0],
            "doi": doi,
            "year": art["Journal"]["JournalIssue"]["PubDate"].get("Year", "")
        })
    return pd.DataFrame(rows)



def _safe_crossref_request(doi: str, retries: int = 3, delay: float = 1.0):
    url = f"https://api.crossref.org/works/{doi}"
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=_XREF_HEADERS, timeout=10)
            if r.status_code == 200:
                # content-type may be text/html for errors, so try json but fall back
                try:
                    return r.json()["message"]
                except (ValueError, KeyError):
                    logging.warning("Crossref non‑JSON body for DOI %s", doi)
                    return None
            elif r.status_code in (404, 410):
                return None  # invalid DOI
            else:           # 429 or 5xx -> retry
                logging.warning("Crossref %s on %s (attempt %d)",
                                r.status_code, doi, attempt+1)
        except requests.RequestException as exc:
            logging.warning("Crossref request error on %s: %s", doi, exc)
        time.sleep(delay)
    return None

def crossref_enrich(df: pd.DataFrame) -> pd.DataFrame:
    """Add licence & PDF link columns; tolerate failures gracefully."""
    rows = []
    for doi in df["doi"].dropna().unique():
        msg = _safe_crossref_request(doi)
        if not msg:
            rows.append({"doi": doi, "license": None, "link_pdf": None})
            continue
        rows.append({
            "doi": doi,
            "license": msg.get("license", [{}])[0].get("URL") if msg.get("license") else None,
            "link_pdf": next((l["URL"] for l in msg.get("link", [])
                              if l.get("content-type") == "application/pdf"), None)
        })
        time.sleep(0.1)  # stay well below polite 50 r/s limit
    return df.merge(pd.DataFrame(rows), on="doi", how="left")
    
def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset="doi").assign(
        sha256=lambda x: x["abstract"].apply(lambda t: hashlib.sha256(t.encode()).hexdigest())
    )

def collect(queries: list[str]) -> pd.DataFrame:
    frames = [pubmed_search(q) for q in queries]
    df = crossref_enrich(pd.concat(frames, ignore_index=True))
    return deduplicate(df)
