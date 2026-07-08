#!/usr/bin/env bash
# Path A — enable CONRAD OpenCL (GPU) on Apple Silicon.
#
# CONRAD's OpenCL fails on arm64 because pyconrad bundles jogamp 2.3.2 (2015),
# whose gluegen has no aarch64 CPU-detection case (throws before any native
# loads). Fix: (1) drop jogamp 2.6.0 (universal arm64 natives) into the pyconrad
# bundle dir so they are globbed onto the classpath before conrad_1.1.0.jar and
# shadow the stale JOCL; (2) build a tiny libOpenCL.dylib that re-exports Apple's
# OpenCL.framework (JOCL's macOS search list doesn't include the framework path);
# conrad_backend.setup() puts .jogamp on DYLD_LIBRARY_PATH before JVM start.
#
# Verified on Apple M1 Max: OpenCLUtil.getStaticContext() -> GPU; CL fan
# projector matches CPU to 0.03% and is ~4000x faster.
set -euo pipefail
cd "$(dirname "$0")/.."

BUNDLE=".venv/lib/python3.12/site-packages/pyconrad/CONRAD 1.1.0"
BASE="https://jogamp.org/deployment/maven"
JARS=(
  "org/jogamp/gluegen/gluegen-rt/2.6.0/gluegen-rt-2.6.0.jar"
  "org/jogamp/gluegen/gluegen-rt/2.6.0/gluegen-rt-2.6.0-natives-macosx-universal.jar"
  "org/jogamp/jocl/jocl/2.6.0/jocl-2.6.0.jar"
  "org/jogamp/jocl/jocl/2.6.0/jocl-2.6.0-natives-macosx-universal.jar"
)

if [ ! -d "$BUNDLE" ]; then echo "pyconrad bundle dir not found: $BUNDLE"; exit 1; fi

echo "Downloading jogamp 2.6.0 (universal arm64 natives) into the pyconrad bundle..."
for j in "${JARS[@]}"; do
  out="$BUNDLE/$(basename "$j")"
  [ -f "$out" ] && { echo "  have $(basename "$j")"; continue; }
  curl -sSL "$BASE/$j" -o "$out"
  echo "  fetched $(basename "$j")"
done

echo "Building OpenCL.framework reexport shim -> .jogamp/libOpenCL.dylib ..."
mkdir -p .jogamp
cat > .jogamp/shim.c <<'EOF'
/* re-export Apple's OpenCL.framework as libOpenCL.dylib for JOCL */
EOF
clang -dynamiclib -arch arm64 -Wl,-reexport_framework,OpenCL -framework OpenCL \
      -o .jogamp/libOpenCL.dylib .jogamp/shim.c
lipo -archs .jogamp/libOpenCL.dylib
echo "Done. Verify:  python -c \"import sys;sys.path.insert(0,'src');import conrad_backend as c;print('OpenCL:',c.opencl_available())\""
