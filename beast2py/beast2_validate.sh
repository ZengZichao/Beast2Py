#!/usr/bin/env bash
# beast2_validate.sh — Headless BEAST2 XML validation.
#
# Validates BEAST2 XML by invoking the BEAST2 v2.7.x XMLParser directly
# (beast.base.parser.XMLParser, via the bundled Beast2Validator tool) without
# requiring JavaFX or the BEAST2 GUI. Prints "VALID: <file>" for each file
# that parses successfully.
#
# Usage:
#   bash beast2_validate.sh <xml-file> [<xml-file> ...]
#
# Environment variables:
#   BEAST2_JAR  path to BEAST.base.jar (default: auto-detected, see below)
#   JAVA_HOME   JDK 17+ home directory. Optional: if unset — or if it points at
#               something older than 17 — a JDK 17+ is searched for among
#               /usr/libexec/java_home, the Homebrew and system JVM directories,
#               and only then `java` on PATH (an old Oracle shim on PATH used to
#               shadow a perfectly usable openjdk@17).
#   BEAST2_VALIDATE_TRACE  set to anything to print which java was chosen
#
# Auto-detection order for BEAST.base.jar:
#   1. $BEAST2_JAR
#   2. newest ~/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar (BEAST2 package dir)
#   3. tools/lib/BEAST.base.jar or lib/BEAST.base.jar next to this script
#
# Exit codes: 0 = all files valid, 1 = at least one file invalid,
#             2 = setup error (jar or Java not found / too old).

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Locate BEAST.base.jar ---------------------------------------------------
if [ -n "${BEAST2_JAR:-}" ] && [ -f "$BEAST2_JAR" ]; then
    :
else
    BEAST2_JAR=""
    for j in "$HOME"/.beast/2.7/BEAST.base/*/lib/BEAST.base.jar; do
        # loop keeps the last (typically newest) match
        [ -f "$j" ] && BEAST2_JAR="$j"
    done
    if [ -z "$BEAST2_JAR" ]; then
        for j in "$SCRIPT_DIR/tools/lib/BEAST.base.jar" "$SCRIPT_DIR/lib/BEAST.base.jar"; do
            if [ -f "$j" ]; then BEAST2_JAR="$j"; break; fi
        done
    fi
fi
if [ -z "${BEAST2_JAR:-}" ] || [ ! -f "$BEAST2_JAR" ]; then
    echo "ERROR: BEAST.base.jar not found." >&2
    echo "  Install BEAST2 v2.7.x, or set BEAST2_JAR=/path/to/BEAST.base.jar" >&2
    exit 2
fi

# --- Locate Java (>= 17, required by BEAST 2.7.8 class files) ----------------
# A JDK 17+ is searched for among *candidates* rather than taken from `$PATH`:
# macOS boxes very often ship an ancient Oracle `java` shim on PATH (1.8.0_501
# reports "java version \"1.8.0_501\"") while a perfectly good openjdk@17 sits
# in /opt/homebrew, unlinked.  Trusting PATH there would fail the third gate on
# a machine that can run BEAST2 just fine, so the first candidate whose major
# version is >= 17 wins and the rest are only reported if none qualifies.
java_major() {
    # echo the major version of a java binary ("8" for 1.8, "17" for 17.0.9)
    local raw
    raw="$("$1" -version 2>&1 | head -1 | sed -E 's/.*version "([^"]*)".*/\1/')"
    case "$raw" in
        1.[0-9]*) echo "$raw" | sed -E 's/^1\.([0-9]+).*/\1/' ;;
        *)        echo "$raw" | sed -E 's/^([0-9]+).*/\1/' ;;
    esac
}

IS_OK_JAVA() {
    local m
    m="$(java_major "$1")"
    [ -n "$m" ] && [ "$m" -ge 17 ] 2>/dev/null
}

CANDIDATES=()
TRACE=""
# BEAST2_JAVA_CANDIDATES: PATH-style list of java binaries to try first, for a
# JDK installed somewhere the probes below do not cover (and for tests).
if [ -n "${BEAST2_JAVA_CANDIDATES:-}" ]; then
    IFS=':' read -r -a _extra_java <<< "$BEAST2_JAVA_CANDIDATES"
    for cand in "${_extra_java[@]}"; do
        [ -x "$cand" ] && CANDIDATES+=("$cand")
    done
fi
[ -n "${JAVA_HOME:-}" ] && [ -x "${JAVA_HOME}/bin/java" ] && CANDIDATES+=("${JAVA_HOME}/bin/java")

if [ -z "${BEAST2_SKIP_JAVA_HOME_TOOL:-}" ] && [ -x /usr/libexec/java_home ]; then
    # macOS: ask for any 17+ JVM (both the "17+" range syntax and the plain "17")
    for spec in "17+" "17"; do
        h="$(/usr/libexec/java_home -v "$spec" 2>/dev/null || true)"
        if [ -n "$h" ] && [ -x "$h/bin/java" ]; then CANDIDATES+=("$h/bin/java"); fi
    done
fi

# Homebrew / Linux / system JVM directories, newest LTS first
# (BEAST2_SKIP_JDK_PROBES=1 disables this block and the glob below, leaving only
#  BEAST2_JAVA_CANDIDATES, JAVA_HOME, /usr/libexec/java_home and PATH — used by
#  the test that checks the "no JDK 17 anywhere" message is actually reachable.)
if [ -z "${BEAST2_SKIP_JDK_PROBES:-}" ]; then
for home in /opt/homebrew/opt/openjdk@21 /opt/homebrew/opt/openjdk@17 \
            /opt/homebrew/opt/openjdk /opt/homebrew/Cellar/openjdk@17 \
            /usr/local/opt/openjdk@21 /usr/local/opt/openjdk@17 /usr/local/opt/openjdk \
            /usr/lib/jvm/java-21-openjdk /usr/lib/jvm/java-17-openjdk \
            /usr/lib/jvm/default-java; do
    for cand in "$home/libexec/openjdk.jdk/Contents/Home/bin/java" "$home/bin/java"; do
        [ -x "$cand" ] && CANDIDATES+=("$cand")
    done
done
for cand in /Library/Java/JavaVirtualMachines/*/Contents/Home/bin/java \
            /usr/lib/jvm/*/bin/java; do
    [ -x "$cand" ] && CANDIDATES+=("$cand")
done
fi

# PATH last: it is the least specific thing available
p="$(command -v java || true)"
[ -n "$p" ] && CANDIDATES+=("$p")

JAVA=""
for cand in "${CANDIDATES[@]}"; do
    if IS_OK_JAVA "$cand"; then JAVA="$cand"; break; fi
    TRACE="$TRACE  rejected (<17 or unreadable): $cand (major $(java_major "$cand"))
"
done

if [ -z "$JAVA" ]; then
    echo "ERROR: BEAST 2.7.8 needs a JDK 17 or newer and none was found." >&2
    printf '%s' "$TRACE" >&2
    echo "  Set JAVA_HOME to a JDK 17+ installation (e.g. brew install openjdk@17)." >&2
    exit 2
fi

JAVA_MAJOR="$(java_major "$JAVA")"
if [ -n "${BEAST2_VALIDATE_TRACE:-}" ]; then
    echo "[beast2_validate] using java $JAVA_MAJOR at $JAVA" >&2
fi

# --- Build classpath ---------------------------------------------------------
# BEAST.base is a self-contained jar; launcher.jar/classes provide the
# beast.pkgmgmt package-manager classes used by Beast2Validator.
CLASSPATH="$BEAST2_JAR:$SCRIPT_DIR/tools/launcher.jar:$SCRIPT_DIR/tools/classes"

# Add-on packages installed under ~/.beast/2.7 (EBSP, BD Skyline Serial, ...),
# skipping the BEAST.base jar already on the classpath.
for j in "$HOME"/.beast/2.7/*/*/lib/*.jar; do
    [ -f "$j" ] || continue
    case "$j" in
        "$BEAST2_JAR") continue ;;
    esac
    case "$j" in
        */BEAST.base.jar) continue ;;
        *) CLASSPATH="$CLASSPATH:$j" ;;
    esac
done

# --- Run the headless validator ----------------------------------------------
if [ -f "$SCRIPT_DIR/tools/classes/Beast2Validator.class" ]; then
    exec "$JAVA" -cp "$CLASSPATH" Beast2Validator "$@"
else
    # Fall back to JDK single-file source mode (compiles in memory)
    exec "$JAVA" -cp "$CLASSPATH" "$SCRIPT_DIR/tools/Beast2Validator.java" "$@"
fi
