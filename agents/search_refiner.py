import openai, json
from config import settings

_SYSTEM = "You are a search strategist for biochemical literature."

_schema = {
    "name": "improved_queries",
    "description": "Return 3‑5 refined queries",
    "parameters": {
        "type": "object",
        "properties": {
            "queries": {
                "type": "array",
                "items": {"type": "string"}
            }
        },
        "required": ["queries"]
    }
}

def improve(initial_queries: list[str]) -> list[str]:
    if not settings.refine_search:
        return initial_queries
    client = openai.OpenAI(api_key=settings.openai_api_key)
    msg = [{"role":"system","content":_SYSTEM},
           {"role":"user","content":str(initial_queries)}]
    res = client.chat.completions.create(
        model=settings.llm_model,
        messages=msg,
        functions=[_schema],
        function_call={"name":"improved_queries"}
    )
    content = res.choices[0].message.function_call.arguments
    return json.loads(content)["queries"]
