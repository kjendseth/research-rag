# agents/collector.py

import logging
import re
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, Medline   # pip install biopython

from config import settings

def _fetch_medline_batch(pmids, retries=3, delay=2.0):
    """
    Fetch a batch of PMIDs in Medline format.
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
                "Medline fetch %s error: %s (retry %d/%d)",
                pmids[:1], e, attempt, retries
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch PMIDs {pmids[:1]} after {retries} retries")

def _fetch_nbib_batch(pmids, retries=3, delay=2.0):
    """
    Fetch the same batch of PMIDs in NBIB format, returning a list of entry strings.
    """
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                rettype="nbib",
                retmode="text"
            )
            raw = handle.read()
            handle.close()
            # Split on lookahead for "PMID- "
            entries = re.split(r'(?=PMID- )', raw)
            # Remove any empty leading chunk
            if entries and not entries[0].startswith("PMID-"):
                entries = entries[1:]
            return entries
        except Exception as e:
            attempt += 1
            logging.warning(
                "NBIB fetch %s error: %s (retry %d/%d)",
                pmids[:1], e, attempt, retries
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed to fetch NBIB for PMIDs {pmids[:1]} after {retries} retries")

def _parse_nbib_entry(entry: str) -> dict[str, list[str]]:
    """
    Parse a single NBIB entry block into a dict of tag -> [values].
    """
    fields: dict[str, list[str]] = {}
    current = None
    for line in entry.splitlines():
        if not line.strip(): 
            continue
        m = re.match(r"^([A-Z0-9]{2,4})- (.*)$", line)
        if m:
            tag, val = m.group(1), m.group(2).strip()
            fields.setdefault(tag, []).append(val)
            current = tag
        elif line.startswith("    ") and current:
            fields[current][-1] += " " + line.strip()
    return fields

def pubmed_search(query: str, cap: int = 5000) -> pd.DataFrame:
    """
    Fetch up to `cap` PubMed records for `query`. Returns a DataFrame
    with top‑level columns plus a full `nbib_fields` dict.
    """
    Entrez.email = settings.ncbi_email

    # 1) total count & cap
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"]); handle.close()
    if total > cap:
        logging.warning("'%s' returned %d hits; capping at %d", query, total, cap)
        total = cap

    # 2) get all PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]; handle.close()

    rows = []
    # 3) process in batches
    batch_size = 500
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]

        # fetch medline + nbib in parallel
        med_records = _fetch_medline_batch(batch)
        nbib_entries = _fetch_nbib_batch(batch)

        if len(nbib_entries) != len(med_records):
            logging.warning(
                "Batch size mismatch: %d Medline vs %d NBIB entries",
                len(med_records), len(nbib_entries)
            )

        for rec, nbib in zip(med_records, nbib_entries):
            pmid        = rec.get("PMID", "")
            nbib_fields = _parse_nbib_entry(nbib)

            # top‑level convenience fields
            title    = " ".join(rec.get("TI","").splitlines())
            abstract = " ".join(rec.get("AB","").splitlines())
            journal  = rec.get("JT","")
            year     = rec.get("DP","")[:4]
            authors  = rec.get("AU",[])
            keywords = rec.get("MH",[])

            # DOI
            doi = next((aid.split()[0] for aid in rec.get("AID",[]) if aid.endswith("[doi]")), None)
            if not doi:
                for val in nbib_fields.get("LID",[]):
                    if "doi" in val.lower():
                        doi = val.split()[0]; break

            # PMC
            pmc = None
            for val in nbib_fields.get("PMC",[]):
                pmc = val.strip(); break

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

        time.sleep(0.5)  # be polite

    return pd.DataFrame(rows)

def collect(queries: list[str]) -> pd.DataFrame:
    """
    Run pubmed_search for each query, concatenate, drop duplicate PMIDs.
    """
    df_list = [pubmed_search(q) for q in queries]
    df = pd.concat(df_list, ignore_index=True)
    return df.drop_duplicates(subset="pmid")