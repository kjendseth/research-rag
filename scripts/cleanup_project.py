#!/usr/bin/env python3
# File: scripts/cleanup_project.py

import sys
from pathlib import Path
import fnmatch

# ─── Ensure project root is on PYTHONPATH ────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.resolve()

import openai
import json
from config import settings

def main():
    import argparse
    p = argparse.ArgumentParser(
        description="Delete vector stores, their attached files, and all local manifest JSONs")
    p.add_argument("project", help="project folder name (e.g. t3cup)")
    args = p.parse_args()

    client   = openai.OpenAI(api_key=settings.openai_api_key)
    project  = args.project
    proj_dir = PROJECT_ROOT / project

    # 1) Attempt to load store IDs (if missing, warn and proceed)
    ids = {}
    sid_file = proj_dir / "store_ids.json"
    if sid_file.exists():
        ids = json.loads(sid_file.read_text())
    else:
        print(f"⚠️ No {sid_file}; skipping vector-store deletion step.")

    # 2) For each store ID, delete attached files then the store
    for vs_id in ids.values():
        print(f"\n→ Cleaning up vector store {vs_id}")
        # 2a) list & delete attached files
        try:
            files = client.vector_stores.files.list(vector_store_id=vs_id).data
        except Exception as e:
            print(f"  ⚠️ Could not list files: {e}")
            files = []
        for f in files:
            fid = getattr(f, "id", None)
            if fid:
                print(f"  • Deleting file {fid}")
                try:
                    client.files.delete(fid)
                except Exception as e:
                    print(f"    ⚠️ {e}")
        # 2b) delete the vector store
        try:
            client.vector_stores.delete(vs_id)
            print(f"  ✓ Deleted vector store {vs_id}")
        except Exception as e:
            print(f"  ⚠️ {e}")

    # 3) Remove ALL local JSON manifests under project folder
    for path in proj_dir.glob("*.json"):
        print(f"✓ Removing {path}")
        path.unlink()

    # 4) Remove root‑level manifest files (.abstract_manifest*.json, .pdf_manifest.json)
    for root_file in PROJECT_ROOT.iterdir():
        if fnmatch.fnmatch(root_file.name, ".abstract_manifest*.json") or root_file.name == ".pdf_manifest.json":
            print(f"✓ Removing {root_file}")
            root_file.unlink()

    print("\n🚀 Cleanup complete.")

if __name__ == "__main__":
    main()