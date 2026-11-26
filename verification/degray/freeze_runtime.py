#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze the local macOS science runtime without overwriting earlier evidence.

FVQual core lookup becomes relative to each archived executable. The SDK and
its dependencies remain external; this is an experiment archive, not a release.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build",type=Path);parser.add_argument("output",type=Path)
    parser.add_argument("--extra-binary",action="append",default=[],help="Relative build path to another science driver")
    args=parser.parse_args();build=args.build.resolve();output=args.output.resolve()
    if platform.system()!="Darwin": raise SystemExit("This local runtime freezer currently supports macOS only")
    if output.exists() and any(output.iterdir()): raise SystemExit("Preserve the existing archive; choose a new output folder")
    binaries=("apps/fvqual-degray","verification/fvqual-circulation-review",*args.extra_binary)
    if any(Path(name).is_absolute() or ".." in Path(name).parts for name in binaries):
        raise SystemExit("Science binary paths must stay within the build")
    for name in (*binaries,"src/libfvqual-core.dylib"):
        target=output/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(build/name,target)
    (output/"src/libfvqual-core.0.dylib").symlink_to("libfvqual-core.dylib")
    for name in binaries:
        subprocess.run(["/usr/bin/install_name_tool","-rpath",str(build/"src"),"@executable_path/../src",str(output/name)],check=True,capture_output=True)
        subprocess.run(["/usr/bin/codesign","--force","--sign","-",str(output/name)],check=True,capture_output=True)
    source=Path(__file__).resolve().parents[2]
    for directory in ("src","include","apps","verification","tests","cmake","python"):
        for path in (source/directory).rglob("*"):
            if path.is_file() and path.suffix in (".cpp",".h",".hpp",".py",".txt",".in") and "output" not in path.parts and "__pycache__" not in path.parts:
                target=output/"sources"/path.relative_to(source);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
    shutil.copy2(source/"CMakeLists.txt",output/"sources/CMakeLists.txt")
    shutil.copytree(build/"generated",output/"generated")
    files={str(path.relative_to(output)):hashlib.sha256(path.read_bytes()).hexdigest()
           for path in output.rglob("*") if path.is_file()}
    metadata={"files":files,"scope":"Frozen FVQual science drivers/core and source; local core lookup; external SDK runtime required",
        "platform":platform.platform(),"original_build":str(build),
        "external_dependencies":{name:subprocess.check_output(["/usr/bin/otool","-L",str(output/name)],text=True) for name in binaries}}
    (output/"archive.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(f"Frozen {len(files)} files in {output}")


if __name__=="__main__": main()
