# agents/collector.py

import requests
import pandas as pd
import time
import hashlib
import json
import logging
import os
from Bio import Entrez  # pip install biopython
from config import settings

# … keep your _safe_crossref_request, crossref_enrich, deduplicate as before …

def pubmed_search(query, retmax=200) -> pd.DataFrame:
    Entrez.email = settings.ncbi_email
    ids = Entrez.read(Entrez.esearch(db="pubmed", term=query,
                                     retmax=retmax))["IdList"]
    fetch = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    articles = Entrez.read(fetch)["PubmedArticle"]
    rows = []

    for art in articles:
        a = art["MedlineCitation"]["Article"]
        # DOI
        doi = next((i for i in a.get("ELocationID", [])
                    if i.attributes["EIdType"] == "doi"), None)
        # Title
        title = a.get("ArticleTitle", "")
        # Abstract
        abstract = " ".join(a.get("Abstract", {}).get("AbstractText", [""]))
        # Journal
        journal = a["Journal"].get("Title", "")
        # Year
        year = a["Journal"]["JournalIssue"]["PubDate"].get("Year", "")
        # Authors
        auths = []
        for au in a.get("AuthorList", []):
            last = au.get("LastName")
            initials = au.get("Initials")
            if last and initials:
                auths.append(f"{last} {initials}")
        # MeSH keywords
        mesh_terms = []
        mesh_list = art["MedlineCitation"].get("MeshHeadingList", [])
        for mh in mesh_list:
            term = mh.get("DescriptorName")
            if term:
                mesh_terms.append(str(term))
        rows.append({
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": auths,
            "keywords": mesh_terms
        })

    return pd.DataFrame(rows)


def collect(queries: list[str]) -> pd.DataFrame:
    frames = [pubmed_search(q) for q in queries]
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates(subset="doi")
    return df