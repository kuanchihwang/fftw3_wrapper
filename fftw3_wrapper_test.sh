#!/usr/bin/env sh

[ -f "${PROGRAM}.py" ] && python3 "${PROGRAM}.py" || exit 1
[ -d "${FFTWPATH}/lib" ] && LD_LIBRARY_PATH="${FFTWPATH}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" "./${PROGRAM}.exe" || exit 2
