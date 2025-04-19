# agents/doi_selector.py
import json, openai
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
    Returns a list of DOIs (strings) proposed for full‑text download.
    """
    client = openai.OpenAI(api_key=settings.openai_api_key)

    resp = client.chat.completions.create(
        model=settings.llm_model,
        messages=[
            {"role": "system", "content": _SYSTEM},
            {"role": "user",   "content": question}
        ],
        tools=[{"type": "file_search"}],
        tool_choice={"type": "file_search"},
        tool_resources={"file_search": {"vector_store_ids": [abstract_store_id], "max_num_results": k}},
        functions=[_FN_SCHEMA],
        function_call={"name": "select_dois"}
    )

    args = json.loads(resp.choices[0].message.function_call.arguments)
    return args["dois"]