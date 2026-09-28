#!/usr/bin/env bash
# Wrapper kept for the repository layout (tests and CI call this path).
# The real script and its Java helper live inside the package so that a
# plain `pip install .` ships a working BEAST2 validator.
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/beast2py/beast2_validate.sh" "$@"
