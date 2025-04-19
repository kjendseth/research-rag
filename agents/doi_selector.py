# agents/doi_selector.py

import json
import openai
from config import settings

_SYSTEM = (
    "You are a literature scout who sees chunks of abstracts and selects the most relevant papers. "
    "Return your answer only as strict JSON in the format:\n"
    '{"dois": ["10.1234/abc", "10.5678/xyz"]}\n'
    "Do not emit any other text."
)

def run(question: str, abstract_store_id: str, k: int = 8) -> list[str]:
    """
    Searches the abstract vector store, then returns a JSON list of DOIs.
    """
    client = openai.OpenAI(api_key=settings.openai_api_key)

    # Retrieval + JSON‑only response via Responses API
    response = client.responses.create(
        model=settings.llm_model,
        input=question,
        instructions=_SYSTEM,
        tools=[{
            "type": "file_search",
            "vector_store_ids": [abstract_store_id],
            "max_num_results": k
        }]
    )

    # Parse the output_text as JSON
    try:
        data = json.loads(response.output_text)
        dois = data.get("dois", [])
        if not isinstance(dois, list):
            raise ValueError("`dois` is not a list")
        return dois
    except (json.JSONDecodeError, ValueError) as e:
        raise RuntimeError(
            f"Failed to parse DOI list from model output: {e}\n"
            f"Raw output: {response.output_text}"
        )