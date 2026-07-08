#!/usr/bin/env python3
"""Create a Zenodo DRAFT deposition for the RabbitCT dataset and upload files.

RUN ON THE LME SERVER (lme31), where the data lives (lme31 has python3 but no
curl/jq). Streams large files (does not load them into RAM). Does NOT publish —
leaves a draft for review. Sandbox dry run: ZENODO_SANDBOX=1 python3 upload_to_zenodo.py

Token: put it in the git-ignored zenodo_secrets.env next to this file (or set
ZENODO_TOKEN in the environment). Never commit real tokens.
"""
import json
import os
import sys
import urllib.request
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = "/disks/data1/share/webdata/fileadmin/Forschung/Software/RabbitCT/download"

# Small files first so any problem surfaces before the multi-GB uploads.
FILES = [
    os.path.join(HERE, "README_data.md"),
    os.path.join(DATA_DIR, "reference_256.vol"),
    os.path.join(DATA_DIR, "rabbitct_develop-v2.zip"),
    os.path.join(DATA_DIR, "rabbitct_512-v2.rctd"),
    os.path.join(DATA_DIR, "rabbitct_1024-v2.rctd"),
]


def load_secrets():
    p = os.path.join(HERE, "zenodo_secrets.env")
    if os.path.isfile(p):
        for line in open(p):
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                os.environ.setdefault(k.strip(), v.strip())


def api(method, url, token, data=None, headers=None, expect=(200, 201)):
    hdrs = {"Authorization": "Bearer " + token}
    if headers:
        hdrs.update(headers)
    req = urllib.request.Request(url, data=data, method=method, headers=hdrs)
    try:
        with urllib.request.urlopen(req) as r:
            body = r.read()
            return r.status, (json.loads(body) if body else {})
    except urllib.error.HTTPError as e:
        return e.code, {"error": e.read().decode("utf-8", "replace")[:400]}


def put_file(bucket, path, token):
    name = os.path.basename(path)
    size = os.path.getsize(path)
    with open(path, "rb") as f:                       # streamed, not read into RAM
        req = urllib.request.Request(bucket + "/" + name, data=f, method="PUT",
                                     headers={"Authorization": "Bearer " + token,
                                              "Content-Type": "application/octet-stream",
                                              "Content-Length": str(size)})
        try:
            with urllib.request.urlopen(req) as r:
                return r.status
        except urllib.error.HTTPError as e:
            print("   ERROR", e.code, e.read().decode("utf-8", "replace")[:300])
            return e.code


def main():
    load_secrets()
    sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
    base = "https://sandbox.zenodo.org/api" if sandbox else "https://zenodo.org/api"
    token = os.environ.get("ZENODO_SANDBOX_TOKEN" if sandbox else "ZENODO_TOKEN") \
        or os.environ.get("ZENODO_TOKEN", "")
    if not token:
        sys.exit("No token. Fill zenodo_secrets.env (see .example) or set ZENODO_TOKEN.")

    print("Endpoint:", base, "(sandbox)" if sandbox else "(PRODUCTION)")
    missing = [f for f in FILES if not os.path.isfile(f)]
    if missing:
        print("WARNING missing files (skipped):"); [print("  ", m) for m in missing]
    files = [f for f in FILES if os.path.isfile(f)]

    print("Creating draft deposition ...")
    st, dep = api("POST", base + "/deposit/depositions", token,
                  data=b"{}", headers={"Content-Type": "application/json"})
    if st not in (200, 201):
        sys.exit("Failed to create deposition: %s %s" % (st, dep))
    dep_id, bucket = dep["id"], dep["links"]["bucket"]
    print("  deposition id:", dep_id)

    for f in files:
        mb = os.path.getsize(f) / 1e6
        print("Uploading %s (%.1f MB) ..." % (os.path.basename(f), mb))
        code = put_file(bucket, f, token)
        print("  -> HTTP", code)
        if code not in (200, 201):
            sys.exit("Upload failed for %s" % f)

    print("Setting metadata ...")
    meta = open(os.path.join(HERE, "zenodo_metadata.json"), "rb").read()
    st, resp = api("PUT", base + "/deposit/depositions/%s" % dep_id, token,
                   data=meta, headers={"Content-Type": "application/json"})
    print("  -> HTTP", st, "" if st == 200 else resp)

    web = base.replace("/api", "")
    print("\nDONE (draft, NOT published).")
    print("Review + publish at: %s/deposit/%s" % (web, dep_id))
    print("To publish after review:")
    print("  python3 - <<'P'\n  import urllib.request\n  urllib.request.urlopen(urllib.request.Request("
          "'%s/deposit/depositions/%s/actions/publish', method='POST', "
          "headers={'Authorization':'Bearer '+open('zenodo_secrets.env').read()}))\n  P" % (base, dep_id))


if __name__ == "__main__":
    main()
