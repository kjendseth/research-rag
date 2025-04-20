# agents/collector.py

import logging
import time
import re
from pathlib import Path

import pandas as pd
from Bio import Entrez, Medline   # pip install biopython

from config import settings

def _fetch_medline_batch(pmids, retries=3, delay=2.0):
    """Fetch Medline‐parsed metadata for a batch of PMIDs."""
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
    raise RuntimeError(f"Failed Medline fetch for {pmids[:1]}")

def _fetch_xml_batch(pmids, retries=3, delay=2.0):
    """
    Fetch the raw PubMed XML for a batch of PMIDs, return list of strings,
    one <PubmedArticle>…</PubmedArticle> block per PMID in order.
    """
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                rettype="xml",
                retmode="text"
            )
            raw_bytes = handle.read()
            handle.close()
            # ensure we have a str
            if isinstance(raw_bytes, (bytes, bytearray)):
                raw = raw_bytes.decode("utf-8")
            else:
                raw = raw_bytes
            # split into article blocks
            parts = re.split(r"(?=<PubmedArticle>)", raw)
            # drop any leading before first tag
            blocks = [p for p in parts if p.startswith("<PubmedArticle>")]
            return blocks
        except Exception as e:
            attempt += 1
            logging.warning(
                "XML fetch %s error: %s (retry %d/%d)",
                pmids[:1], e, attempt, retries
            )
            time.sleep(delay)
    raise RuntimeError(f"Failed XML fetch for {pmids[:1]}")

def pubmed_search(query: str, cap: int = 5000) -> pd.DataFrame:
    """
    Fetch up to `cap` PubMed records for `query`. Returns a DataFrame with:
      pmid, doi, pmc, title, abstract, journal, year, authors, keywords, xml
    """
    Entrez.email = settings.ncbi_email

    # 1) total count & cap
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"])
    handle.close()
    if total > cap:
        logging.warning("'%s' returned %d hits; capping at %d", query, total, cap)
        total = cap

    # 2) retrieve PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]
    handle.close()

    rows = []
    batch_size = 500
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        med_meta = _fetch_medline_batch(batch)
        xml_blocks = _fetch_xml_batch(batch)

        if len(xml_blocks) != len(med_meta):
            logging.warning(
                "Batch size mismatch: %d Medline vs %d XML blocks",
                len(med_meta), len(xml_blocks)
            )

        for rec, xml in zip(med_meta, xml_blocks):
            pmid    = rec.get("PMID", "")
            title   = " ".join(rec.get("TI","").splitlines())
            abstract= " ".join(rec.get("AB","").splitlines())
            journal = rec.get("JT","")
            year    = rec.get("DP","")[:4]
            authors = rec.get("AU",[])
            keywords= rec.get("MH",[])

            # DOI
            doi = next((a.split()[0] for a in rec.get("AID",[]) if a.endswith("[doi]")), None)

            # PMC from XML
            pmc = None
            m = re.search(r'<ArticleId IdType="pmc">(PMC\d+)</ArticleId>', xml)
            if m:
                pmc = m.group(1)

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
                "xml": xml
            })

        time.sleep(0.5)  # be polite to NCBI

    return pd.DataFrame(rows)

def collect(queries: list[str]) -> pd.DataFrame:
    """
    Run pubmed_search for each query term, concatenate results,
    and drop duplicate PMIDs.
    """
    df_list = [pubmed_search(q) for q in queries]
    df = pd.concat(df_list, ignore_index=True)
    return df.drop_duplicates(subset="pmid")