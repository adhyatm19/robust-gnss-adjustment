#!/usr/bin/env bash
# Build a pinned RTKLIB rnx2rtkp in an isolated, gitignored location.
#
# Produces two binaries from the SAME pinned source tree:
#   tools/rtklib/bin/rnx2rtkp-stock  unmodified RTKLIB 2.4.3 b34
#   tools/rtklib/bin/rnx2rtkp-hp     identical except for the output-precision
#                                    patch third_party/rtklib/rnx2rtkp-output-precision.patch
#                                    (printf widths only; no numerical change)
# and tools/rtklib/build_info.json recording commit, patch hash and compiler.
#
# Nothing is installed system-wide. Requires git, gcc, make.
# Usage: scripts/build_rtklib.sh [TOOLS_DIR]   (default: <repo>/tools)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="${1:-$REPO_ROOT/tools}"
RTK_DIR="$TOOLS_DIR/rtklib"
RTKLIB_URL="https://github.com/tomojitakasu/RTKLIB"
RTKLIB_TAG="v2.4.3-b34"
RTKLIB_COMMIT="180043ee24b6d2b168f98b64be15f69d50046b1a"
PATCH="$REPO_ROOT/third_party/rtklib/rnx2rtkp-output-precision.patch"

mkdir -p "$RTK_DIR/bin"
if [ ! -d "$RTK_DIR/src/.git" ]; then
    git clone --quiet --depth 1 --branch "$RTKLIB_TAG" "$RTKLIB_URL" "$RTK_DIR/src"
fi
actual="$(git -C "$RTK_DIR/src" rev-parse HEAD)"
if [ "$actual" != "$RTKLIB_COMMIT" ]; then
    echo "ERROR: RTKLIB checkout is $actual, expected pinned $RTKLIB_COMMIT" >&2
    exit 1
fi
if [ -n "$(git -C "$RTK_DIR/src" status --porcelain --untracked-files=no)" ]; then
    echo "ERROR: pinned RTKLIB source tree has local modifications" >&2
    exit 1
fi

build_variant() {  # $1 = variant name, $2 = patch file or empty
    local work="$RTK_DIR/build-$1"
    rm -rf "$work"
    mkdir -p "$work"
    git -C "$RTK_DIR/src" archive HEAD src app/consapp/rnx2rtkp | tar -x -C "$work"
    if [ -n "$2" ]; then
        patch --quiet -p1 -d "$work" < "$2"
    fi
    # The b34 makefile links lib/iers (only needed with -DIERS_MODEL, which the
    # default OPTS do not set); override LDLIBS instead of editing the makefile.
    make --quiet -C "$work/app/consapp/rnx2rtkp/gcc" LDLIBS="-lm -lrt" >/dev/null
    cp "$work/app/consapp/rnx2rtkp/gcc/rnx2rtkp" "$RTK_DIR/bin/rnx2rtkp-$1"
}

build_variant stock ""
build_variant hp "$PATCH"

cat > "$RTK_DIR/build_info.json" <<EOF
{
  "rtklib_url": "$RTKLIB_URL",
  "rtklib_tag": "$RTKLIB_TAG",
  "rtklib_commit": "$RTKLIB_COMMIT",
  "patch_file": "third_party/rtklib/rnx2rtkp-output-precision.patch",
  "patch_sha256": "$(sha256sum "$PATCH" | cut -d' ' -f1)",
  "make_overrides": "LDLIBS=-lm -lrt",
  "compiler": "$(gcc --version | head -1)",
  "stock_sha256": "$(sha256sum "$RTK_DIR/bin/rnx2rtkp-stock" | cut -d' ' -f1)",
  "hp_sha256": "$(sha256sum "$RTK_DIR/bin/rnx2rtkp-hp" | cut -d' ' -f1)"
}
EOF
echo "built: $RTK_DIR/bin/rnx2rtkp-stock and $RTK_DIR/bin/rnx2rtkp-hp"
