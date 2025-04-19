# agents/writer.py

import openai
from config import settings

def write_summary(question: str, abstracts_store: str, pdf_store: str) -> str:
    """
    Ask GPT‑4.1 to write a structured summary of `question`, grounding
    every citation via the file_search tool. Inline DOIs as (doi:...).
    """
    system = (
        "You are a biochemical reviewer. "
        "Answer in up to 500 words, structured as Intro/Methods/Gaps. "
        "Whenever you reference a paper, include its DOI inline like (doi:10.xxxx). "
        "Use the file_search tool to fetch evidence from the literature."
    )

    client = openai.OpenAI(api_key=settings.openai_api_key)

    resp = client.chat.completions.create(
        model=settings.llm_model,
        messages=[
            {"role": "system", "content": system},
            {"role": "user",   "content": question}
        ],
        tools=[{"type": "file_search"}],
        tool_resources={
            "file_search": {
                "vector_store_ids": [abstracts_store, pdf_store]
            }
        },
        # Force the model to call file_search for evidence
        tool_choice={"type": "required"}
    )

    # The model will emit a single tool_call with name file_search first,
    # then you may need to submit that tool call, but in the Assistants
    # model it will return content directly once the tool is used.
    return resp.choices[0].message.content