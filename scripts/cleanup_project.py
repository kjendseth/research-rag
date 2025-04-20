#!/usr/bin/env python3
# File: scripts/cleanup_project.py

import sys
import os
from pathlib import Path
import fnmatch
import json
import openai

# ─── Ensure project root is on PYTHONPATH ────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT))

def main():
    import argparse
    p = argparse.ArgumentParser(
        description="Delete vector stores, their attached files, and all local manifest JSONs")
    p.add_argument("project", help="project folder name (e.g. t3cup)")
    args = p.parse_args()

    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("Error: OPENAI_API_KEY not set in environment.")
        sys.exit(1)

    client = openai.OpenAI(api_key=api_key)
    proj_dir = PROJECT_ROOT / args.project

    # 1) Load store IDs if present
    sid_file = proj_dir / "store_ids.json"
    ids = {}
    if sid_file.exists():
        ids = json.loads(sid_file.read_text())
    else:
        print(f"⚠️ No {sid_file}; skipping vector‐store deletion step.")

    # 2) Delete each vector store and its attached files
    for vs_id in ids.values():
        print(f"\n→ Cleaning up vector store {vs_id}")

        # delete attached files
        try:
            attached = client.vector_stores.files.list(vector_store_id=vs_id).data
        except Exception as e:
            print(f"  ⚠️ Could not list files: {e}")
            attached = []
        for f in attached:
            fid = getattr(f, "id", None)
            if fid:
                print(f"  • Deleting file {fid}")
                try:
                    client.files.delete(fid)
                except Exception as e:
                    print(f"    ⚠️ {e}")

        # delete the vector store
        try:
            client.vector_stores.delete(vs_id)
            print(f"  ✓ Deleted vector store {vs_id}")
        except Exception as e:
            print(f"  ⚠️ {e}")

    # 3) Remove all JSON manifests in the project folder
    if proj_dir.exists():
        for path in proj_dir.glob("*.json"):
            print(f"✓ Removing {path}")
            path.unlink()

    # 4) Remove root‐level manifest files (.abstract_manifest*.json, .pdf_manifest.json)
    for root_file in PROJECT_ROOT.iterdir():
        if fnmatch.fnmatch(root_file.name, ".abstract_manifest*.json") or root_file.name == ".pdf_manifest.json":
            print(f"✓ Removing {root_file}")
            root_file.unlink()

    print("\n🚀 Cleanup complete.")

if __name__ == "__main__":
    main()