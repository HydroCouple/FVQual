#!/bin/bash
# S3.4 falsifier — the fvqual run subcommand.
set -u
cd "$(dirname "$0")/.."
SRC=apps/fvqual/runcommand.cpp
BAK=$SRC.s34bak
BUILD="cmake --build build -j 8"
RUN="./build/tests/fvqual_tests --gtest_filter=FVQualRunCli.*"

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

# M1 the manifest is never written
mutate 's/if \(!recorder\.finishToFile\(manifestPath, message,\n                               IO::RunStatus::Completed\)\)/if (false)/s' \
    && verdict 1 || echo "M1: NOT APPLIED"

# M2 snapshots are captured but never delivered
mutate 's/if \(captured && !ioThread\.enqueue\(std::move\(snapshot\)\)\)/if (false)/' \
    && verdict 2 || echo "M2: NOT APPLIED"

# M3 instances lose the document id
mutate 's/std::make_unique<\n                                             FVQual::Component::FVQualComponent>\(\n                                             id, &info\)/info.createComponentInstance()/s' \
    && verdict 3 || echo "M3: NOT APPLIED"

# M4 the document's writers are ignored
mutate 's/for \(const IO::WriterSpec &writer : spec\.writers\)/for (const IO::WriterSpec \&writer : std::vector<IO::WriterSpec>{})/' \
    && verdict 4 || echo "M4: NOT APPLIED"
