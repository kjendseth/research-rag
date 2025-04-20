#!/usr/bin/env python3
"""
scripts/delete_all_project_files.py

Deletes all uploaded files in your OpenAI account that look like
JSON/JSONL or PDF uploads for this project, looping until none remain.
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

def main():
    print("Starting deletion of JSON/JSONL/PDF files…")
    while True:
        # 2) Fetch one page of files (up to 100)
        resp = client.files.list(limit=100)
        to_delete = []
        for f in resp.data:
            # filename attribute may vary by SDK version
            name = getattr(f, "filename", None) or getattr(f, "name", None) or f.id
            if name.endswith((".json", ".jsonl", ".pdf")):
                to_delete.append((f.id, name))

        if not to_delete:
            print("No more matching files to delete.")
            break

        print(f"Found {len(to_delete)} files to delete on this pass:")
        for fid, fname in to_delete:
            try:
                client.files.delete(fid)
                print(f"✓ Deleted {fname} ({fid})")
            except Exception as e:
                print(f"✗ Failed to delete {fname} ({fid}): {e}")

    print("\nAll matching files have been removed.")

if __name__ == "__main__":
    main()