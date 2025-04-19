# main.py  (replace entire file)

import argparse, time
from config import get_store_ids
from agents import collector, indexer, downloader, doi_selector

def pubmed_ingest(search_term: str, project: str):
    """
    1️⃣  Harvest abstracts via PubMed and push them into the project’s
        abstracts vector‑store.  No PDFs yet.
    """
    ids = get_store_ids(project)
    df  = collector.collect([search_term])
    indexer.index_abstracts(df, ids["abstracts"])
    print("✅ PubMed ingestion done.")

def enrich_fulltext(question: str, project: str):
    """
    2️⃣  Run DOI‑selection over the abstracts store, download PDFs, embed them.
    """
    ids  = get_store_ids(project)
    dois = doi_selector.run(question, ids["abstracts"])
    print(f"🔍 Selected {len(dois)} DOIs.")
    for doi in dois:
        pdf = downloader.fetch_pdf(doi, dest_dir=f"{project}/pdfs")
        if pdf:
            indexer.index_pdf(pdf, ids["pdfs"])
        time.sleep(1)   # be polite
    print("🎉 Full‑text enrichment done.")

if __name__ == "__main__":
    p = argparse.ArgumentParser(prog="research‑rag")
    sub = p.add_subparsers(dest="cmd", required=True)

    ing = sub.add_parser("ingest", help="PubMed → abstracts store")
    ing.add_argument("search_term")
    ing.add_argument("project")

    enr = sub.add_parser("enrich", help="Select DOIs and fetch PDFs")
    enr.add_argument("project")
    enr.add_argument("question")

    a = p.parse_args()
    if a.cmd == "ingest":
        pubmed_ingest(a.search_term, a.project)
    elif a.cmd == "enrich":
        enrich_fulltext(a.question, a.project)