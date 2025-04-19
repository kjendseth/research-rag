# agents/collector.py

import requests
import pandas as pd
import time
import hashlib
import json
import logging
import os
from Bio import Entrez          # pip install biopython
from config import settings

# … keep your _safe_crossref_request and crossref_enrich if you want enrichment …

def pubmed_search(query, cap=5000) -> pd.DataFrame:
    """
    Fetch all records for `query`, up to `cap`. Extract pmid, doi,
    title, abstract, journal, year, authors, keywords.
    """
    Entrez.email = settings.ncbi_email

    # 1) get total count
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"])
    if total > cap:
        print(f"⚠️ {total} hits for '{query}', capping to {cap}")
        total = cap

    # 2) fetch all PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]

    # 3) fetch in batches
    records = []
    batch = 500
    for i in range(0, len(pmids), batch):
        chunk = pmids[i : i + batch]
        recs = _fetch_medline_batch(chunk)
        records.extend(recs)
        time.sleep(0.5)

    # 4) parse out fields
    rows = []
    for art in records:
        pmid = art.get("PMID", "")
        title = art.get("TI", "").replace("\n", " ")
        abstract = art.get("AB", "").replace("\n", " ")
        journal = art.get("JT", "")
        year = art.get("DP", "")[:4]
        # authors
        auths = art.get("AU", [])
        # keywords (MeSH)
        mesh = art.get("MH", [])

        # DOI from AID field
        doi = next(
            (aid.split()[0] for aid in art.get("AID", [])
             if aid.endswith("[doi]")),
            None
        )

        rows.append({
            "pmid": pmid,
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": auths,
            "keywords": mesh
        })

    return pd.DataFrame(rows)


def _fetch_medline_batch(pmids, retries=3, delay=2.0):
    """
    Helper to efetch and parse Medline text, with retries on IncompleteRead.
    """
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                rettype="medline",
                retmode="text"
            )
            records = list(Entrez.parse(handle))
            handle.close()
            return records
        except Exception as e:
            attempt += 1
            logging.warning("Fetch batch %s error: %s (retry %d/%d)",
                            pmids[:1], e, attempt, retries)
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch PMIDs {pmids[:1]} after {retries} retries")


def collect(queries: list[str]) -> pd.DataFrame:
    """
    Run pubmed_search over each query and concatenate unique PMIDs.
    """
    dfs = [pubmed_search(q) for q in queries]
    df = pd.concat(dfs, ignore_index=True)

    # dedupe by pmid (so we keep all, even if doi is missing)
    df = df.drop_duplicates(subset="pmid")
    return df