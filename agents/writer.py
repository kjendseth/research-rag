import openai, json
from config import settings

def write_summary(question: str, abstracts_store: str, pdf_store: str) -> str:
    system = (
        "You are a biochemical reviewer.\n"
        "Answer the question in ≤500 words.\n"
        "When you refer to a paper, embed its DOI inline like (doi:10.1016/j.jmb.2025.01.015).\n"
        "Always ground statements with file_search from the provided vector stores."
    )
    tools = [
        {"type": "file_search", "vector_store_id": abstracts_store},
        {"type": "file_search", "vector_store_id": pdf_store}
    ]
    client = openai.OpenAI(api_key=settings.openai_api_key)
    res = client.chat.completions.create(
        model=settings.llm_model,
        messages=[{"role":"system","content":system},
                  {"role":"user","content":question}],
        tools=tools
    )
    return res.choices[0].message.content
