import requests, pandas as pd, time, hashlib, json, openai
from Bio import Entrez          # pip install biopython
from config import settings

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

def crossref_enrich(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for d in df["doi"].dropna().unique():
        r = requests.get(f"https://api.crossref.org/works/{d}").json()["message"]
        rows.append({
            "doi": d,
            "license": r.get("license", [{}])[0].get("URL"),
            "link_pdf": next((l["URL"] for l in r.get("link", [])
                              if l["content-type"]=="application/pdf"), None)
        })
        time.sleep(0.1)
    return df.merge(pd.DataFrame(rows), on="doi", how="left")

def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset="doi").assign(
        sha256=lambda x: x["abstract"].apply(lambda t: hashlib.sha256(t.encode()).hexdigest())
    )

def collect(queries: list[str]) -> pd.DataFrame:
    frames = [pubmed_search(q) for q in queries]
    df = crossref_enrich(pd.concat(frames, ignore_index=True))
    return deduplicate(df)
