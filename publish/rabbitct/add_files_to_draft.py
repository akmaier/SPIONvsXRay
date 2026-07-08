#!/usr/bin/env python3
"""Add more files to an existing Zenodo DRAFT deposition (run on lme31).

Used to append the ORIGINAL (non-v2) RabbitCT files to the draft that already
holds the -v2 versions, so the record has both. Does NOT publish.

Usage:  python3 add_files_to_draft.py <deposition_id>
Token from zenodo_secrets.env (git-ignored) or ZENODO_TOKEN.
"""
import json
import os
import sys
import urllib.request
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = "/disks/data1/share/webdata/fileadmin/Forschung/Software/RabbitCT/download"

# The updated README (mentions both versions; overwrites the one on the draft)
# plus the ORIGINAL (2011) versions to add alongside the already-uploaded -v2 files.
EXTRA_FILES = [
    os.path.join(HERE, "README_data.md"),
    os.path.join(DATA_DIR, "rabbitct_develop.zip"),
    os.path.join(DATA_DIR, "rabbitct_512.rctd"),
    os.path.join(DATA_DIR, "rabbitct_1024.rctd"),
]


def load_secrets():
    p = os.path.join(HERE, "zenodo_secrets.env")
    if os.path.isfile(p):
        for line in open(p):
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                os.environ.setdefault(k.strip(), v.strip())


def put_file(bucket, path, token):
    name = os.path.basename(path)
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        req = urllib.request.Request(bucket + "/" + name, data=f, method="PUT",
                                     headers={"Authorization": "Bearer " + token,
                                              "Content-Type": "application/octet-stream",
                                              "Content-Length": str(size)})
        try:
            with urllib.request.urlopen(req) as r:
                return r.status
        except urllib.error.HTTPError as e:
            print("   ERROR", e.code, e.read().decode("utf-8", "replace")[:300]); return e.code


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: add_files_to_draft.py <deposition_id>")
    dep_id = sys.argv[1]
    load_secrets()
    base = "https://sandbox.zenodo.org/api" if os.environ.get("ZENODO_SANDBOX") == "1" else "https://zenodo.org/api"
    token = os.environ.get("ZENODO_TOKEN", "")
    if not token:
        sys.exit("No token (zenodo_secrets.env / ZENODO_TOKEN).")

    req = urllib.request.Request(base + "/deposit/depositions/" + dep_id,
                                 headers={"Authorization": "Bearer " + token})
    with urllib.request.urlopen(req) as r:
        dep = json.loads(r.read())
    bucket = dep["links"]["bucket"]
    print("Adding to deposition", dep_id)
    for f in EXTRA_FILES:
        if not os.path.isfile(f):
            print("  MISSING (skip):", f); continue
        print("Uploading %s (%.1f MB) ..." % (os.path.basename(f), os.path.getsize(f) / 1e6))
        print("  -> HTTP", put_file(bucket, f, token))
    print("Done (still a draft, NOT published). Review at",
          base.replace("/api", "") + "/deposit/" + dep_id)


if __name__ == "__main__":
    main()
