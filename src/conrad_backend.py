"""CONRAD/pyconrad backend bootstrap.

Handles the Apple-Silicon environment quirks discovered at M0:
  * CONRAD needs Java 8; JPype >= 1.6 dropped Java 8 -> we pin JPype1 1.5.0.
  * On arm64 the JVM must also be arm64 -> we ship a local arm64 Zulu JDK 8
    under .jdk8/ (see scripts/install_jdk8.sh) and point JAVA_HOME at it.

Import this and call ``setup()`` once before using any CONRAD class.
"""
from __future__ import annotations
import glob
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
_JDK_GLOB = str(REPO_ROOT / ".jdk8" / "zulu*" / "Contents" / "Home")

_started = False


def _find_local_jdk8() -> str | None:
    """Return JAVA_HOME of the bundled arm64 Zulu JDK 8, if present."""
    hits = sorted(glob.glob(_JDK_GLOB))
    for h in hits:
        if (Path(h) / "jre" / "lib" / "jli" / "libjli.dylib").exists() or \
           (Path(h) / "lib" / "jli" / "libjli.dylib").exists() or \
           (Path(h) / "bin" / "java").exists():
            return h
    return None


def ensure_java_home() -> str:
    """Point JAVA_HOME at a suitable JDK 8 (prefer the bundled arm64 one)."""
    local = _find_local_jdk8()
    if local:
        os.environ["JAVA_HOME"] = local
        return local
    jh = os.environ.get("JAVA_HOME")
    if jh:
        return jh
    raise RuntimeError(
        "No JDK 8 found. Run scripts/install_jdk8.sh (Apple Silicon) or set "
        "JAVA_HOME to a Java 8 install."
    )


def setup(max_ram: str = "8G"):
    """Start the JVM and initialise CONRAD (idempotent). Returns the pyconrad module."""
    global _started
    ensure_java_home()
    import pyconrad  # imported after JAVA_HOME is set
    if not _started:
        pyconrad.setup_pyconrad(max_ram=max_ram)
        _started = True
    return pyconrad


def class_getter(package: str):
    """Convenience: pyconrad.ClassGetter for a CONRAD package (after setup)."""
    import pyconrad
    return pyconrad.ClassGetter(package)


if __name__ == "__main__":
    jh = ensure_java_home()
    print("JAVA_HOME:", jh)
    pc = setup()
    print("pyconrad", getattr(pc, "__version__", "?"), "JVM up.")
    cg = class_getter("edu.stanford.rsl.conrad.data.numeric")
    g = cg.Grid2D(4, 4)
    g.setAtIndex(1, 1, 2.5)
    ok = g.getAtIndex(1, 1) == 2.5
    print("CONRAD bridge OK:", ok)
    sys.exit(0 if ok else 1)
