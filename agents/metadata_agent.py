# agents/metadata_agent.py

import json
import openai
from config import settings
from agents.metadata_tools import semantic_search, search_by_author, search_by_pmc

FUNCTIONS = [
    {
        "name": "semantic_search",
        "description": "Find papers by semantic search over abstracts/title",
        "parameters": {
            "type": "object",
            "properties": {
                "vector_store_id": {"type": "string"},
                "query": {"type": "string"},
                "limit": {"type": "integer", "default": 100}
            },
            "required": ["vector_store_id", "query"]
        }
    },
    {
        "name": "search_by_author",
        "description": "Find all papers co‑authored by someone",
        "parameters": {
            "type": "object",
            "properties": {
                "vector_store_id": {"type": "string"},
                "author_substr": {"type": "string"},
                "limit": {"type": "integer", "default": 1000}
            },
            "required": ["vector_store_id", "author_substr"]
        }
    },
    {
        "name": "search_by_pmc",
        "description": "Find all papers with a PMC ID",
        "parameters": {
            "type": "object",
            "properties": {
                "vector_store_id": {"type": "string"},
                "limit": {"type": "integer", "default": 1000}
            },
            "required": ["vector_store_id"]
        }
    }
]

def agent_query(question: str, store_id: str) -> str:
    client = openai.OpenAI(api_key=settings.openai_api_key)

    # 1) Ask the model which function to call (it won't know your store_id)
    resp = client.chat.completions.create(
        model=settings.llm_model,
        messages=[
            {"role": "system", "content": "You can call functions to fetch papers by metadata or semantic search."},
            {"role": "user",   "content": question}
        ],
        functions=FUNCTIONS,
        function_call="auto"
    )

    msg = resp.choices[0].message

    # If no function was chosen, just return the text
    if not msg.function_call:
        return msg.content

    # 2) Parse the function call and inject your real store_id
    fn_name = msg.function_call.name
    args = json.loads(msg.function_call.arguments)
    args["vector_store_id"] = store_id

    # 3) Execute the correct tool
    if fn_name == "semantic_search":
        results = semantic_search(**args)
    elif fn_name == "search_by_author":
        results = search_by_author(**args)
    elif fn_name == "search_by_pmc":
        results = search_by_pmc(**args)
    else:
        return f"Error: function {fn_name} not implemented."

    # 4) Hand the raw results list back to GPT for a human‑friendly summary
    follow = client.chat.completions.create(
        model=settings.llm_model,
        messages=[
            {"role": "system",    "content": "Here are the function results you requested."},
            {"role": "assistant", "content": json.dumps(results)},
            {"role": "user",      "content": "Please summarize these entries, listing each DOI and a one‑sentence snippet."}
        ]
    )
    return follow.choices[0].message.content