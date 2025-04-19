# agents/writer.py

import openai
from config import settings

def write_summary(question: str, store_id: str) -> str:
    """
    Uses the Responses API to run a file_search over the given vector store
    and returns a grounded summary with inline DOIs.
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
            "vector_store_ids": [store_id],
            "max_num_results": 8
        }]
    )

    return response.output_text