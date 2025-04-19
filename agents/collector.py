# agents/collector.py

import logging
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, Medline  # pip install biopython

from config import settings

def _fetch_medline_batch(pmids, retries=3, delay=2.0):
    """
    Helper to efetch and parse a batch of PMIDs from PubMed in Medline format,
    retrying on errors. Returns a list of Medline record dicts.
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
            records = list(Medline.parse(handle))
            handle.close()
            return records
        except Exception as e:
            attempt += 1
            logging.warning(
                "Fetch batch %s error: %s (retry %d/%d)",
                pmids[:1], e, attempt, retries
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch PMIDs {pmids[:1]} after {retries} retries")


def pubmed_search(query: str, cap: int = 5000) -> pd.DataFrame:
    """
    Fetch up to `cap` PubMed records for `query`. Parses out:
    - pmid
    - doi (if present)
    - title, abstract, journal, year
    - authors list
    - MeSH keywords list
    Returns a pandas DataFrame.
    """
    Entrez.email = settings.ncbi_email

    # 1) Get total count
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"])
    handle.close()
    if total > cap:
        logging.warning("Query '%s' returned %d hits; capping at %d", query, total, cap)
        total = cap

    # 2) Retrieve all PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]
    handle.close()

    # 3) Fetch records in batches
    records = []
    batch_size = 500
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        recs = _fetch_medline_batch(batch)
        records.extend(recs)
        time.sleep(0.5)  # be polite to NCBI

    # 4) Parse fields
    rows = []
    for rec in records:
        pmid     = rec.get("PMID", "")
        title    = rec.get("TI", "").replace("\n", " ")
        abstract = rec.get("AB", "").replace("\n", " ")
        journal  = rec.get("JT", "")
        year     = rec.get("DP", "")[:4]
        authors  = rec.get("AU", [])               # list of strings
        keywords = rec.get("MH", [])               # MeSH headings list

        # Extract DOI if present
        doi = next(
            (aid.split()[0] for aid in rec.get("AID", []) if aid.endswith("[doi]")),
            None
        )

        rows.append({
            "pmid": pmid,
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": authors,
            "keywords": keywords
        })

    return pd.DataFrame(rows)


def collect(queries: list[str]) -> pd.DataFrame:
    """
    Run pubmed_search for each query in `queries`, concatenate the results,
    and drop duplicate PMIDs. Returns a unified DataFrame.
    """
    df_list = []
    for q in queries:
        df_q = pubmed_search(q)
        df_list.append(df_q)
    df = pd.concat(df_list, ignore_index=True)
    # Drop duplicate records by PMID, keep first occurrence
    df = df.drop_duplicates(subset="pmid")
    return df