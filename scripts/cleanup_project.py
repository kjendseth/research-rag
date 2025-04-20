#!/usr/bin/env python3
# File: scripts/cleanup_project.py

import sys
from pathlib import Path

# ─── Ensure project root is on PYTHONPATH ────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT))

import openai
import json
from config import settings

def main():
    import argparse
    p = argparse.ArgumentParser(
        description="Delete vector stores, their attached files, and all local manifest JSONs for a project")
    p.add_argument("project", help="project folder name (e.g. t3cup)")
    args = p.parse_args()

    client   = openai.OpenAI(api_key=settings.openai_api_key)
    project  = args.project
    proj_dir = PROJECT_ROOT / project
    sid_file = proj_dir / "store_ids.json"

    # 1) Load vector‑store IDs if present
    ids = {}
    if sid_file.exists():
        ids = json.loads(sid_file.read_text())
    else:
        print(f"⚠️ Warning: no {sid_file}; skipping vector‐store deletion.")

    # 2) For each store, delete all its attached files then the store itself
    for store_key, vs_id in ids.items():
        print(f"\n→ Cleaning up {store_key} store: {vs_id}")
        # 2a) list & delete attached files
        try:
            attached = client.vector_stores.files.list(vector_store_id=vs_id).data
        except Exception as e:
            print(f"  ⚠️ Could not list files for {vs_id}: {e}")
            attached = []
        for att in attached:
            fid = getattr(att, "file_id", None) or getattr(att, "id", None)
            if not fid:
                continue
            print(f"  • Deleting file {fid}")
            try:
                client.files.delete(fid)
            except Exception as e:
                print(f"    ⚠️ Error deleting file {fid}: {e}")
        # 2b) delete the vector store
        try:
            client.vector_stores.delete(vs_id)
            print(f"  ✓ Deleted vector store {vs_id}")
        except Exception as e:
            print(f"  ⚠️ Error deleting vector store {vs_id}: {e}")

    # 3) Remove any local JSON manifests
    patterns = [
        proj_dir / "store_ids.json",
        PROJECT_ROOT / ".abstract_manifest_*.json",
        PROJECT_ROOT / ".abstract_manifest.json",
        PROJECT_ROOT / ".pdf_manifest.json"
    ]
    # also catch any stray project‑specific JSON
    patterns.append(PROJECT_ROOT / f"{project}" / "store_ids.json")

    for pat in patterns:
        for path in PROJECT_ROOT.glob(pat.name) if "*" in pat.name else [pat]:
            if path.exists():
                path.unlink()
                print(f"✓ Removed local {path}")

    print("\n🚀 Cleanup complete.")

if __name__ == "__main__":
    main()