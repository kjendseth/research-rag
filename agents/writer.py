# agents/writer.py

import openai
from config import settings

def write_summary(question: str, pdf_store: str) -> str:
    """
    Uses the Responses API to run a file_search over the full‑text vector store
    and then returns a grounded 500‑word summary with inline DOIs.
    """
    system = (
        "You are a biochemical reviewer. Answer in up to 500 words, "
        "structured with headings (Introduction, Methods, Gaps). "
        "When citing papers, include inline DOIs like (doi:10.xxxx)."
    )

    client = openai.OpenAI(api_key=settings.openai_api_key)

    response = client.responses.create(
        model=settings.llm_model,
        input=question,
        instructions=system,
        tools=[{
            "type": "file_search",
            "vector_store_ids": [pdf_store],    # <-- only full‑text store
            "max_num_results": 8
        }]
    )

    return response.output_text