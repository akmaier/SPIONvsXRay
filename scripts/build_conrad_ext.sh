#!/usr/bin/env bash
# Compile our CONRAD extension (fixed OpenCL fan backprojector with pixel spacing)
# against the pyconrad CONRAD jar. Output classes + the .cl resource go to
# conrad_ext/out/, which conrad_backend puts on the classpath (CONRAD_DEV_DIRS)
# BEFORE the jar so the fixed class is available. See DEVLOG (CL backprojector fix).
set -euo pipefail
cd "$(dirname "$0")/.."

JH="$(ls -d .jdk8/zulu*/Contents/Home 2>/dev/null | head -1)"
[ -n "$JH" ] || { echo "arm64 JDK 8 not found (run scripts/install_jdk8.sh)"; exit 1; }
JAR="$(ls ".venv/lib/python3.12/site-packages/pyconrad/CONRAD 1.1.0/conrad_1.1.0.jar")"

SRC="conrad_ext/edu/stanford/rsl/tutorial/fan"
OUT="conrad_ext/out"
rm -rf "$OUT"; mkdir -p "$OUT"
# Ensure the sibling kernel the class also loads (Ray) is present next to Pixel.
CSRC=/Users/maier/Documents/CONRAD/src/edu/stanford/rsl/tutorial/fan
[ -f "$SRC/FanBeamBackProjectorRay.cl" ] || cp "$CSRC/FanBeamBackProjectorRay.cl" "$SRC/" 2>/dev/null || true
echo "Compiling patched FanBeamBackprojector2D (spacing fix) against CONRAD jar ..."
"$JH/bin/javac" -source 8 -target 8 -cp "$JAR" -d "$OUT" "$SRC/FanBeamBackprojector2D.java"
# .cl resources must sit next to the .class for getResourceAsStream()
cp "$SRC"/*.cl "$OUT/edu/stanford/rsl/tutorial/fan/"
echo "Built:"; find "$OUT" -type f
