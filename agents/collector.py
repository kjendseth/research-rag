# agents/collector.py

import logging
import re
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, Medline   # pip install biopython

from config import settings

def _fetch_medline_batch(pmids, retries=3, delay=2.0):
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
                "Fetch Medline batch %s error: %s (retry %d/%d)",
                pmids[:1], e, attempt, retries
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch PMIDs {pmids[:1]} after {retries} retries")

def _parse_nbib_entry(entry: str) -> dict[str, list[str]]:
    """
    Parse a single NBIB entry block into a dict of tag -> [values].
    Handles multi-line continuations indented by spaces.
    """
    fields: dict[str, list[str]] = {}
    current_tag = None
    for line in entry.splitlines():
        if not line.strip():
            continue
        m = re.match(r"^([A-Z0-9]{2,4})- (.*)$", line)
        if m:
            tag, val = m.group(1), m.group(2).strip()
            fields.setdefault(tag, []).append(val)
            current_tag = tag
        elif line.startswith("    ") and current_tag:
            fields[current_tag][-1] += " " + line.strip()
    return fields

def pubmed_search(query: str, cap: int = 5000) -> pd.DataFrame:
    """
    Fetch up to `cap` PubMed records for `query`. Returns a DataFrame with:
      - top‑level columns: pmid, doi, pmc, title, abstract, journal, year, authors, keywords
      - nbib_fields: dict of *all* raw NBIB tags for each record
    """
    Entrez.email = settings.ncbi_email

    # 1) total count & cap
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"]); handle.close()
    if total > cap:
        logging.warning("Query '%s' returned %d hits; capping at %d", query, total, cap)
        total = cap

    # 2) get PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]; handle.close()

    # 3) fetch Medline metadata in batches
    med_records = []
    for i in range(0, len(pmids), 500):
        batch = pmids[i : i + 500]
        med_records.extend(_fetch_medline_batch(batch))
        time.sleep(0.5)

    # 4) fetch all NBIB text
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="nbib", retmode="text")
    raw_nbib = handle.read().strip(); handle.close()
    nbib_entries = raw_nbib.split("\n\n")

    # 5) assemble rows
    rows = []
    for rec, nbib in zip(med_records, nbib_entries):
        pmid        = rec.get("PMID", "")
        nbib_fields = _parse_nbib_entry(nbib)

        # convenience fields
        title    = " ".join(rec.get("TI", "").splitlines())
        abstract = " ".join(rec.get("AB", "").splitlines())
        journal  = rec.get("JT", "")
        year     = rec.get("DP", "")[:4]
        authors  = rec.get("AU", [])
        keywords = rec.get("MH", [])

        # DOI extraction
        doi = next((aid.split()[0] for aid in rec.get("AID", []) if aid.endswith("[doi]")), None)
        if not doi:
            for val in nbib_fields.get("LID", []):
                if "doi" in val.lower():
                    doi = val.split()[0]
                    break

        # PMC extraction
        pmc = None
        for val in nbib_fields.get("PMC", []):
            # entries may appear as "PMC1234567"
            pmc = val.strip()
            break

        rows.append({
            "pmid": pmid,
            "doi": doi,
            "pmc": pmc,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": authors,
            "keywords": keywords,
            "nbib_fields": nbib_fields
        })

    return pd.DataFrame(rows)

def collect(queries: list[str]) -> pd.DataFrame:
    """
    Run pubmed_search for each query term, concatenate results,
    and drop duplicate PMIDs.
    """
    dfs = [pubmed_search(q) for q in queries]
    df = pd.concat(dfs, ignore_index=True)
    return df.drop_duplicates(subset="pmid")