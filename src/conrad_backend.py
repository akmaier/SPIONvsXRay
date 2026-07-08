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
OPENCL_AVAILABLE = None   # set True/False after setup() probes it


def _enable_opencl_env():
    """Point dyld at the local OpenCL reexport shim (Path A GPU enablement).

    The jogamp 2.6.0 jars (universal arm64 natives) live in the pyconrad bundle
    dir so they're globbed onto the classpath before conrad_1.1.0.jar, shadowing
    the stale 2.3.2 JOCL. Apple's OpenCL.framework isn't on JOCL's macOS search
    list, so .jogamp/libOpenCL.dylib re-exports it; dyld finds it via
    DYLD_LIBRARY_PATH. Must be set before the JVM starts. See scripts/install_opencl.sh.
    """
    shim_dir = REPO_ROOT / ".jogamp"
    if (shim_dir / "libOpenCL.dylib").exists():
        cur = os.environ.get("DYLD_LIBRARY_PATH", "")
        parts = [p for p in cur.split(":") if p]
        if str(shim_dir) not in parts:
            os.environ["DYLD_LIBRARY_PATH"] = ":".join([str(shim_dir)] + parts)
        return True
    return False


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


def _enable_conrad_ext():
    """Put our compiled CONRAD extension (fixed CL fan backprojector) on the
    classpath BEFORE conrad_1.1.0.jar via CONRAD_DEV_DIRS (pyconrad prepends it).
    Built by scripts/build_conrad_ext.sh."""
    out = REPO_ROOT / "conrad_ext" / "out"
    if (out / "edu" / "stanford" / "rsl" / "tutorial" / "fan").exists():
        cur = os.environ.get("CONRAD_DEV_DIRS", "")
        parts = [p for p in cur.split(";") if p]
        if str(out) not in parts:
            os.environ["CONRAD_DEV_DIRS"] = ";".join([str(out)] + parts)
        return True
    return False


def setup(max_ram: str = "8G"):
    """Start the JVM and initialise CONRAD (idempotent). Returns the pyconrad module."""
    global _started
    ensure_java_home()
    _enable_opencl_env()                 # before JVM start
    _enable_conrad_ext()                 # before JVM start
    import pyconrad  # imported after JAVA_HOME is set
    if not _started:
        pyconrad.setup_pyconrad(max_ram=max_ram)
        _started = True
    return pyconrad


def opencl_available() -> bool:
    """Probe (once) whether CONRAD OpenCL initialises on this machine (Path A)."""
    global OPENCL_AVAILABLE
    if OPENCL_AVAILABLE is None:
        setup()
        try:
            OCU = class_getter("edu.stanford.rsl.conrad.opencl").OpenCLUtil
            OCU.getStaticContext()
            OPENCL_AVAILABLE = True
        except Exception:
            OPENCL_AVAILABLE = False
    return OPENCL_AVAILABLE


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
