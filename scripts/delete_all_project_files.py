#!/usr/bin/env python3
"""
scripts/delete_all_project_files.py

Deletes all uploaded files in your OpenAI account that look like
JSON/JSONL or PDF uploads for this project.
"""
import os
import sys
import openai

# 1) Initialize
api_key = os.getenv("OPENAI_API_KEY")
if not api_key:
    print("Error: set OPENAI_API_KEY")
    sys.exit(1)
client = openai.OpenAI(api_key=api_key)

# 2) Fetch one page of files (up to 100)
print("Fetching file list…")
resp = client.files.list(limit=100)
to_delete = []
for f in resp.data:
    name = getattr(f, "filename", None) or f.id
    if name.endswith((".json", ".jsonl", ".pdf")):
        to_delete.append((f.id, name))

# 3) Delete them
print(f"Found {len(to_delete)} files to delete.")
for fid, fname in to_delete:
    try:
        client.files.delete(fid)
        print(f"✓ Deleted {fname} ({fid})")
    except Exception as e:
        print(f"✗ Failed to delete {fname} ({fid}): {e}")

print("\nAll matching files on this page have been removed.")