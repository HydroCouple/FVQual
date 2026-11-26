#!/bin/bash
# S3.3 falsifier — FVQual's met series and kinetics deck arguments.
set -u
cd "$(dirname "$0")/.."
SRC=src/component/fvqualcomponent.cpp
BAK=$SRC.s33bak
BUILD="cmake --build build -j 8"
RUN="./build/tests/fvqual_tests --gtest_filter=FVQualForcingArguments.*"

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

# M1 the met argument never enables meteorology
mutate 's/        !m_solver\.meteorologyEnabled\(\)\)\n        m_solver\.enableMeteorology\(\);/        false)\n        m_solver.enableMeteorology();/s' \
    && verdict 1 || echo "M1: NOT APPLIED"

# M2 the met series is never applied in update
mutate 's/m_solver\.meteorology\(\)\.airTemperature\.setUniform\(air\);/(void)air;/' \
    && verdict 2 || echo "M2: NOT APPLIED"

# M3 the series is read at the wrong clock (seconds as days)
mutate 's/m_solver\.time\(\) \/ 86400\.0\);/m_solver.time());/' \
    && verdict 3 || echo "M3: NOT APPLIED"

# M4 a kinetics deck parses but registers nothing
mutate 's/if \(m_hasKinetics && !m_reactionFile\.registerSpecies\(m_registry, message\)\)/if (false)/' \
    && verdict 4 || echo "M4: NOT APPLIED"

# M5 rate laws never configured into the solver
mutate 's/if \(!m_reactionFile\.configure\(m_solver, message\)\)/if (false)/' \
    && verdict 5 || echo "M5: NOT APPLIED"

# M6 a malformed deck is swallowed
mutate 's/if \(!m_reactionFile\.read\(kineticsPath, message\)\)\n            return false;/(void)m_reactionFile.read(kineticsPath, message);/s' \
    && verdict 6 || echo "M6: NOT APPLIED"
