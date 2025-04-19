# Research‑RAG 🧬📚

A modular Retrieval‑Augmented Generation (RAG) pipeline that
automates **biochemical literature discovery, indexing and analysis**.

> **Tech stack**  
> OpenAI GPT‑4/4o · OpenAI Vector Stores · PubMed/WoS/Crossref · Unpaywall  
> Playwright (optional EZproxy PDFs) · pandas·pypdf·Streamlit · UMAP/Plotly

---

## Features

| Stage | What it does |
|-------|--------------|
| **Search Refiner** (toggle) | GPT‑4.1 improves your raw keyword list into precise PubMed/WoS queries. |
| **Collector** | Fetches titles, abstracts & DOIs from PubMed, Web of Science, Crossref. |
| **Downloader** | Grabs PDFs via Unpaywall, Crossref full‑text links or headless EZproxy. |
| **Indexer** | Creates **two** OpenAI vector stores: abstracts & full text (chunked). |
| **Analyst** | Answers questions with `file_search` across both stores—cites DOIs only. |
| **Writer** | Produces structured summaries / proposal boiler‑plates with inline citations. |
| **Dashboard** | Streamlit UI: search, UMAP embedding map, chunk inspector, export BibTeX. |

---

## Installation

git clone https://github.com/kjendseth/research-rag.git
cd research-rag
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m playwright install chromium   # one‑time, optional

export OPENAI_API_KEY=sk‑...
export NCBI_EMAIL=you@example.com
# optional:
export WOS_API_KEY=...
export EZPROXY_PREFIX=https://ezproxy.youruni.no/login?url=

---

## Quick Start

# 1. Build / refresh vector stores
python main.py --project atpase --keywords "P‑type ATPase mechanism"

# 2. Ask a research question
python main.py --project atpase --question "What open questions remain about E2P formation?"

# 3. Interactive exploration
streamlit run dashboards/app.py


## Running after installation

cd path/to/research-rag
source .venv/bin/activate          # or .\.venv\Scripts\Activate.ps1 on Windows

---

## Testing

pytest -q

---

## Directory Layout

```text
research‑rag/
├── agents/          ← modular agent logic
│   ├── search_refiner.py
│   ├── collector.py
│   ├── downloader.py
│   ├── indexer.py
│   ├── analyst.py
│   └── writer.py
├── dashboards/
│   └── app.py       ← Streamlit UI
├── tests/
│   └── test_pipeline.py
├── config.py        ← central settings / feature flags
├── main.py          ← CLI orchestrator
├── requirements.txt
└── README.md


