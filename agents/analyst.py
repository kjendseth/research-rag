import openai, json
from config import settings

def search(question: str, vstore_id: str) -> list[dict]:
    tool_call = {
        "type": "file_search",
        "vector_store_id": vstore_id,
        "query": question,
        "k": 8
    }
    res = openai.OpenAI(api_key=settings.openai_api_key).chat.completions.create(
        model=settings.llm_model,
        messages=[{"role":"user","content":question}],
        tools=[tool_call],
        tool_choice={"type":"file_search"}
    )
    return res.choices[0].message.tool_calls[0].function.arguments
