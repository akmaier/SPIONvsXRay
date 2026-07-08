#!/usr/bin/env bash
# Rebuild a patched conrad_1.1.0.jar from the CONRAD source checkout by recompiling
# only the changed classes against the released jar and splicing them (+ updated .cl
# resources) back in with `jar uf`. Avoids a full 1205-file Eclipse build.
#
# Prereqs: arm64 Zulu JDK 8 in ./.jdk8 (scripts/install_jdk8.sh), pyconrad installed
# (provides the base jar), and a CONRAD source checkout with the fixes committed.
# Output: publish/conrad/conrad_1.1.0.jar  and  publish/conrad/CONRAD_1.1.0.zip
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(pwd)"
CONRAD="${CONRAD_SRC:-$HOME/Documents/CONRAD}"
[ -d "$CONRAD/src" ] || { echo "CONRAD source not found at $CONRAD (set CONRAD_SRC)"; exit 1; }

JH="$(ls -d .jdk8/zulu*/Contents/Home 2>/dev/null | head -1)"
[ -n "$JH" ] || { echo "arm64 JDK 8 not found (run scripts/install_jdk8.sh)"; exit 1; }
BASEJAR="$(ls .venv/lib/python3.12/site-packages/pyconrad/'CONRAD 1.1.0'/conrad_1.1.0.jar)"

# The fixed source files (compiled) and resource files (copied as-is). All
# compiled .class outputs (incl. nested classes) are spliced automatically.
JAVA=(
  "edu/stanford/rsl/conrad/data/numeric/NumericGridOperator.java"
  "edu/stanford/rsl/conrad/data/numeric/NumericPointwiseOperators.java"
  "edu/stanford/rsl/conrad/data/numeric/opencl/OpenCLGridOperators.java"
  "edu/stanford/rsl/tutorial/fan/FanBeamBackprojector2D.java"
)
RES=(
  "edu/stanford/rsl/conrad/data/numeric/opencl/PointwiseOperators.cl"
  "edu/stanford/rsl/tutorial/fan/FanBeamBackProjectorPixel.cl"
)

OUT="$(mktemp -d)"; trap 'rm -rf "$OUT"' EXIT
CP="$BASEJAR:$(printf "$CONRAD/lib/%s:" $(cd "$CONRAD/lib" && ls *.jar))"
echo "Compiling changed classes against released jar ..."
"$JH/bin/javac" -source 8 -target 8 -encoding UTF-8 -cp "$CP" -d "$OUT" \
  $(for f in "${JAVA[@]}"; do echo "$CONRAD/src/$f"; done)
for r in "${RES[@]}"; do mkdir -p "$OUT/$(dirname "$r")"; cp "$CONRAD/src/$r" "$OUT/$r"; done

mkdir -p publish/conrad
cp "$BASEJAR" publish/conrad/conrad_1.1.0.jar
# splice every compiled class + every resource, relative to $OUT
( cd "$OUT" && "$JH/bin/jar" uf "$REPO/publish/conrad/conrad_1.1.0.jar" \
    $(find . -name '*.class' | sed 's#^\./##') \
    $(for r in "${RES[@]}"; do echo "$r"; done) )
echo "Patched jar: publish/conrad/conrad_1.1.0.jar (sha1 $(shasum -a1 publish/conrad/conrad_1.1.0.jar | cut -c1-12))"

# Optional: rebuild the FAU release zip if the canonical one is present alongside.
if [ -f publish/conrad/CONRAD_1.1.0.orig.zip ]; then
  echo "Repackaging CONRAD_1.1.0.zip from canonical zip ..."
  S="$(mktemp -d)"; ( cd "$S" && unzip -q "$REPO/publish/conrad/CONRAD_1.1.0.orig.zip" -x "__MACOSX/*" \
    && cp "$REPO/publish/conrad/conrad_1.1.0.jar" "CONRAD 1.1.0/conrad_1.1.0.jar" \
    && zip -q -r -X "$REPO/publish/conrad/CONRAD_1.1.0.zip" "CONRAD 1.1.0" )
  rm -rf "$S"; echo "Release zip: publish/conrad/CONRAD_1.1.0.zip"
fi
echo "Done. See publish/conrad/RELEASE.md for hosting steps."
