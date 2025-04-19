# agents/doi_selector.py

import json
import openai
from config import settings

_SYSTEM = (
    "You are a literature scout working only with abstracts.\n"
    "Analyse the retrieved chunks and return up to 20 DOIs in JSON:\n"
    '{"dois": ["10.1234/abc", "10.5678/xyz"]}'
)

_FN_SCHEMA = {
    "name": "select_dois",
    "parameters": {
        "type": "object",
        "properties": {
            "dois": {
                "type": "array",
                "items": {"type": "string"}
            }
        },
        "required": ["dois"]
    }
}

def run(question: str, abstract_store_id: str, k: int = 8) -> list[str]:
    """
    Searches the abstract vector store, then calls the select_dois function
    to get a JSON list of DOIs.
    """
    client = openai.OpenAI(api_key=settings.openai_api_key)

    # 1) Retrieval + function call via Responses API
    response = client.responses.create(
        model=settings.llm_model,
        input=question,
        instructions=_SYSTEM,
        tools=[{
            "type": "file_search",
            "vector_store_id": abstract_store_id,
            "max_num_results": k
        }],
        functions=[_FN_SCHEMA],
        function_call={"name": "select_dois"}
    )

    # 2) Extract the arguments from the tool call
    # Depending on SDK version, the tool call may live in .tool_calls or .function_call
    if hasattr(response, "tool_calls") and response.tool_calls:
        func = response.tool_calls[0].function
        args = json.loads(func.arguments)
    elif hasattr(response, "function_call"):
        args = json.loads(response.function_call.arguments)
    else:
        raise RuntimeError("No function call returned by DOI selector")

    return args["dois"]