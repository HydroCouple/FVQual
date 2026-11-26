#!/bin/bash
# S3.2 falsifier — FVQual's mesh argument.
set -u
cd "$(dirname "$0")/.."
SRC=src/component/fvqualcomponent.cpp
BAK=$SRC.s32bak
BUILD="cmake --build build -j 8"
RUN="./build/tests/fvqual_tests --gtest_filter=FVQualMeshArgument.*:FVQualPluginTest.*"

cp "$SRC" "$BAK"
restore() { cp "$BAK" "$SRC"; }
trap 'restore; rm -f "$BAK"; touch "$SRC"; $BUILD >/dev/null 2>&1' EXIT

verdict() {
    if ! $BUILD >/dev/null 2>&1; then echo "M$1: BUILD FAILED"; return; fi
    if $RUN >/dev/null 2>&1; then echo "M$1: SURVIVED"; else echo "M$1: CAUGHT"; fi
}

mutate() {
    restore
    perl -0pi -e "$1" "$SRC"
    if cmp -s "$SRC" "$BAK"; then return 1; fi
    touch "$SRC"
    return 0
}

# M1 a supplied mesh no longer sizes the columns
mutate 's/if \(m_meshArgument != nullptr && m_meshArgument->hasMesh\(\)\)\n        m_columnCount =\n            static_cast<int>\(m_meshArgument->mesh\(\)\.faceCount\(\)\);/\/\/ mutant: columns stay the rectangle'"'"'s/s' \
    && verdict 1 || echo "M1: NOT APPLIED"

# M2 prepare always builds the rectangle
mutate 's/const bool suppliedMesh =\n        m_meshArgument != nullptr && m_meshArgument->hasMesh\(\);/const bool suppliedMesh = false;/s' \
    && verdict 2 || echo "M2: NOT APPLIED"

# M3 a supplied mesh's bed is flattened to -depth anyway
mutate 's/if \(suppliedMesh && !definition\.nodeZ\.empty\(\)\)/if (false)/' \
    && verdict 3 || echo "M3: NOT APPLIED"

# M4 the default rectangle degenerates
mutate 's/rectangularBasin\(m_columnCount, m_length, m_width\)/rectangularBasin(1, m_length, m_width)/' \
    && verdict 4 || echo "M4: NOT APPLIED"

# M5 the per-face mean divides by the wrong count
mutate 's/sum \/ static_cast<double>\(end - begin\);/sum \/ static_cast<double>(end - begin + 1);/' \
    && verdict 5 || echo "M5: NOT APPLIED"
