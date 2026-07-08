#!/usr/bin/env bash
# Create a Zenodo DRAFT deposition for the RabbitCT dataset and upload the files.
# RUN THIS ON THE LME SERVER (lme31), where the data lives and bandwidth is good.
#
# It does NOT publish — it leaves a draft for you to review + click Publish on
# Zenodo. Test on sandbox first: ZENODO_SANDBOX=1 ./upload_to_zenodo.sh
#
# Requirements: bash, curl, jq. Token (do NOT paste in chat):
#   export ZENODO_TOKEN=<your Zenodo personal access token: deposit:write,deposit:actions>
set -euo pipefail

: "${ZENODO_TOKEN:?Set ZENODO_TOKEN (Zenodo personal access token) first}"
BASE="https://zenodo.org/api"
[ "${ZENODO_SANDBOX:-0}" = "1" ] && BASE="https://sandbox.zenodo.org/api"

HERE="$(cd "$(dirname "$0")" && pwd)"
META="$HERE/zenodo_metadata.json"
README="$HERE/README_data.md"

# RabbitCT data on the LME server (edit if paths differ):
DATA_DIR="/disks/data1/share/webdata/fileadmin/Forschung/Software/RabbitCT/download"
FILES=(
  "$DATA_DIR/rabbitct_512.rctd"
  "$DATA_DIR/rabbitct_1024.rctd"
  "$DATA_DIR/reference_256.vol"
  "$DATA_DIR/rabbitct_develop.zip"
  "$README"
)

auth=(-H "Authorization: Bearer $ZENODO_TOKEN")

echo "Zenodo endpoint: $BASE"
echo "Creating draft deposition ..."
resp="$(curl -sS "${auth[@]}" -H "Content-Type: application/json" -X POST "$BASE/deposit/depositions" -d '{}')"
dep_id="$(echo "$resp" | jq -r '.id')"
bucket="$(echo "$resp" | jq -r '.links.bucket')"
[ "$dep_id" != "null" ] || { echo "Failed to create deposition:"; echo "$resp"; exit 1; }
echo "  deposition id: $dep_id"
echo "  bucket: $bucket"

for f in "${FILES[@]}"; do
  [ -f "$f" ] || { echo "  MISSING: $f (skipping)"; continue; }
  name="$(basename "$f")"
  sz="$(du -h "$f" | cut -f1)"
  echo "Uploading $name ($sz) ..."
  code="$(curl -sS "${auth[@]}" --upload-file "$f" "$bucket/$name" -o /dev/null -w '%{http_code}')"
  echo "  -> HTTP $code"
done

echo "Setting metadata ..."
code="$(curl -sS "${auth[@]}" -H "Content-Type: application/json" -X PUT \
  "$BASE/deposit/depositions/$dep_id" -d @"$META" -o /tmp/zenodo_meta_resp.json -w '%{http_code}')"
echo "  -> HTTP $code"

echo
echo "DONE (draft, NOT published)."
echo "Review + publish at: ${BASE%/api}/deposit/$dep_id"
echo "When happy, publish with:"
echo "  curl -H \"Authorization: Bearer \$ZENODO_TOKEN\" -X POST $BASE/deposit/depositions/$dep_id/actions/publish"
