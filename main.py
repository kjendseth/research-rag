import argparse
from config import get_store_ids
from agents import search_refiner, collector, downloader, indexer, writer

def pipeline(project, keywords):
    ids = get_store_ids(project)
    queries = search_refiner.improve(keywords)
    df = collector.collect(queries)
    indexer.index_abstracts(df, ids["abstracts"])
    for doi in df["doi"].dropna():
        pdf = downloader.fetch_pdf(doi, dest_dir=f"{project}/pdfs")
        if pdf:
            indexer.index_pdf(pdf, ids["pdfs"])

def rag_cli(project, question):
    ids = get_store_ids(project)
    answer = writer.write_summary(question, ids["pdfs"])
    print(answer)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--project", required=True,
                   help="folder / namespace for this literature corpus")
    p.add_argument("--keywords", nargs="+",
                   help="search keywords (run pipeline)")
    p.add_argument("--question",
                   help="ask a question over existing stores")
    a = p.parse_args()

    if a.keywords:
        pipeline(a.project, a.keywords)
    if a.question:
        rag_cli(a.project, a.question)
