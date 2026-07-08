#!/usr/bin/env bash
# Install a native arm64 (Apple Silicon) Azul Zulu JDK 8 into ./.jdk8/
# (git-ignored, no sudo). CONRAD needs Java 8, and on arm64 the JVM must be
# arm64 to be loadable by an arm64 Python via JPype. See DEVLOG (M0).
set -euo pipefail
cd "$(dirname "$0")/.."

DEST=".jdk8"
URL="https://cdn.azul.com/zulu/bin/zulu8.94.0.17-ca-fx-jdk8.0.492-macosx_aarch64.tar.gz"

if ls "$DEST"/zulu*/Contents/Home/bin/java >/dev/null 2>&1; then
  echo "arm64 JDK 8 already present under $DEST/"; exit 0
fi

mkdir -p "$DEST"
echo "Downloading Zulu JDK 8 (arm64) ..."
curl -sL "$URL" -o "$DEST/zulu8.tar.gz"
tar xzf "$DEST/zulu8.tar.gz" -C "$DEST"
rm -f "$DEST/zulu8.tar.gz"
JH="$(ls -d "$DEST"/zulu*/Contents/Home)"
echo "Installed JDK 8 at: $JH"
"$JH/bin/java" -version
