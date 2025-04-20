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
        description="Delete vector stores, all their attached files, and local manifests for a project")
    p.add_argument("project", help="project folder name (e.g. t3cup)")
    args = p.parse_args()

    client = openai.OpenAI(api_key=settings.openai_api_key)
    project = args.project
    proj_dir = PROJECT_ROOT / project
    sid_file = proj_dir / "store_ids.json"

    if not sid_file.exists():
        print(f"Error: no {sid_file}; nothing to clean up.")
        return

    # 1) Load vector‑store IDs
    ids = json.loads(sid_file.read_text())

    for store_key, vs_id in ids.items():
        print(f"\n→ Cleaning up {store_key} store: {vs_id}")

        # 2) List and delete attached files
        try:
            attached = client.vector_stores.files.list(vector_store_id=vs_id).data
        except Exception as e:
            print(f"  ⚠️ Could not list files for {vs_id}: {e}")
            attached = []

        for att in attached:
            # VectorStoreFile.id is the file-upload ID
            fid = getattr(att, "file_id", None) or getattr(att, "id", None)
            if not fid:
                print(f"  ⚠️ Could not determine file ID for attachment: {att}")
                continue
            print(f"  • Deleting file {fid}")
            try:
                client.files.delete(fid)
            except Exception as e:
                print(f"    ⚠️ Error deleting file {fid}: {e}")

        # 3) Delete the vector store itself
        try:
            client.vector_stores.delete(vs_id)
            print(f"  ✓ Deleted vector store {vs_id}")
        except Exception as e:
            print(f"  ⚠️ Error deleting vector store {vs_id}: {e}")

    # 4) Remove local store_ids.json
    sid_file.unlink()
    print(f"\n✓ Removed local {sid_file}")

    # 5) Remove local manifest files
    abstract_manifest = PROJECT_ROOT / f".abstract_manifest_{ids['abstracts']}.json"
    pdf_manifest      = PROJECT_ROOT / ".pdf_manifest.json"
    for m in (abstract_manifest, pdf_manifest):
        if m.exists():
            m.unlink()
            print(f"✓ Removed local manifest {m}")

    print("\n🚀 Cleanup complete.")

if __name__ == "__main__":
    main()