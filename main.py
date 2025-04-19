# main.py

import argparse, time
from config import get_store_ids
from agents import collector, indexer, downloader, doi_selector, writer

def pubmed_ingest(search_term: str, project: str):
    ids = get_store_ids(project)
    df  = collector.collect([search_term])
    indexer.index_abstracts(df, ids["abstracts"])
    print("✅ PubMed ingestion done.")

def enrich_fulltext(project: str, question: str):
    ids  = get_store_ids(project)
    dois = doi_selector.run(question, ids["abstracts"])
    print(f"🔍 Selected {len(dois)} DOIs:")
    for doi in dois:
        print("  •", doi)
    for doi in dois:
        pdf = downloader.fetch_pdf(doi, dest_dir=f"{project}/pdfs")
        if pdf:
            indexer.index_pdf(pdf, ids["pdfs"])
        time.sleep(1)
    print("🎉 Full‑text enrichment done.")

def query_store(project: str, question: str, store: str):
    """
    Run RAG over either 'abstracts' or 'pdfs' store, printing a summary.
    """
    ids = get_store_ids(project)
    if store not in ids:
        raise ValueError(f"Unknown store '{store}'. Must be one of {list(ids.keys())}.")
    summary = writer.write_summary(question, ids[store])
    print("\n" + summary)

if __name__ == "__main__":
    p = argparse.ArgumentParser(prog="research‑rag")
    sub = p.add_subparsers(dest="cmd", required=True)

    ing = sub.add_parser("ingest", help="PubMed → abstracts store")
    ing.add_argument("search_term")
    ing.add_argument("project")

    enr = sub.add_parser("enrich", help="Select DOIs & harvest PDFs")
    enr.add_argument("project")
    enr.add_argument("question")

    qry = sub.add_parser("query", help="Run RAG over abstracts or PDFs")
    qry.add_argument("project")
    qry.add_argument("question")
    qry.add_argument(
        "--store", choices=["abstracts", "pdfs"], default="pdfs",
        help="Which vector store to query (default: pdfs)"
    )

    args = p.parse_args()

    if args.cmd == "ingest":
        pubmed_ingest(args.search_term, args.project)
    elif args.cmd == "enrich":
        enrich_fulltext(args.project, args.question)
    elif args.cmd == "query":
        query_store(args.project, args.question, args.store)