# agents/collector.py

import logging
import time
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
                db="pubmed", id=",".join(pmids),
                rettype="medline", retmode="text"
            )
            records = list(Medline.parse(handle))
            handle.close()
            return records
        except Exception as e:
            attempt += 1
            logging.warning("Medline fetch %s error: %s (retry %d/%d)",
                            pmids[:1], e, attempt, retries)
            time.sleep(delay)
    raise RuntimeError(f"Failed Medline fetch for {pmids[:1]}")

def _fetch_xml_batch(pmids, retries=3, delay=2.0):
    """Fetch the raw PubMed XML for a batch of PMIDs, return list of strings."""
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(
                db="pubmed", id=",".join(pmids),
                rettype="xml", retmode="text"
            )
            full = handle.read()
            handle.close()
            # split out each <PubmedArticle>…</PubmedArticle> block
            parts = full.split("</PubmedArticle>")
            blocks = [p + "</PubmedArticle>" for p in parts if "<PubmedArticle>" in p]
            return blocks
        except Exception as e:
            attempt += 1
            logging.warning("XML fetch %s error: %s (retry %d/%d)",
                            pmids[:1], e, attempt, retries)
            time.sleep(delay)
    raise RuntimeError(f"Failed XML fetch for {pmids[:1]}")

def pubmed_search(query: str, cap: int = 5000) -> pd.DataFrame:
    """
    Fetch up to `cap` hits for `query`. Returns DataFrame with columns:
      pmid, doi, pmc, title, abstract, journal, year, authors, keywords, xml
    """
    Entrez.email = settings.ncbi_email

    # 1) count & cap
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    total = int(Entrez.read(handle)["Count"]); handle.close()
    if total > cap:
        logging.warning("'%s' → %d hits (capping to %d)", query, total, cap)
        total = cap

    # 2) collect PMIDs
    handle = Entrez.esearch(db="pubmed", term=query, retmax=total)
    pmids = Entrez.read(handle)["IdList"]; handle.close()

    rows = []
    batch = 500
    for i in range(0, len(pmids), batch):
        chunk = pmids[i : i + batch]
        med_meta = _fetch_medline_batch(chunk)
        xml_blocks = _fetch_xml_batch(chunk)

        if len(xml_blocks) != len(med_meta):
            logging.warning(
                "Batch size mismatch: %d Medline vs %d XML blocks",
                len(med_meta), len(xml_blocks)
            )

        for rec, xml in zip(med_meta, xml_blocks):
            pmid    = rec.get("PMID","")
            title   = " ".join(rec.get("TI","").splitlines())
            abstract= " ".join(rec.get("AB","").splitlines())
            journal = rec.get("JT","")
            year    = rec.get("DP","")[:4]
            authors = rec.get("AU",[])
            keywords= rec.get("MH",[])

            # DOI
            doi = next((a.split()[0] for a in rec.get("AID",[]) if a.endswith("[doi]")), None)
            # PMC
            pmc = None
            for block in [xml]:
                m = re.search(r"<ArticleId IdType=\"pmc\">(PMC\d+)</ArticleId>", block)
                if m:
                    pmc = m.group(1)
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
                "xml": xml
            })

        time.sleep(0.5)

    return pd.DataFrame(rows)

def collect(queries: list[str]) -> pd.DataFrame:
    dfs = [pubmed_search(q) for q in queries]
    df = pd.concat(dfs, ignore_index=True)
    return df.drop_duplicates(subset="pmid")