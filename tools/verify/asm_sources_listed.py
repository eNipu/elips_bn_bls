#!/usr/bin/env python3
"""
Every .S in src/arith/ must appear in every build that links the library.

WHY THIS EXISTS. There are three source lists, and they drift. CMakeLists.txt
builds the C library, bindings/python/build_ffi.py builds the Python
extension, and they are maintained by hand. Adding fp_addsub_x86_64.S to the
first and not the second produced a wheel that imported and then died with

    ImportError: undefined symbol: elips_fp_sub_6_x86_64

which CI caught, but only after the commit, and only because the Python job
happens to import the module. A file added for a curve or a platform nobody
tests on that path would not have been caught at all.

The wasm build (bindings/js/build.sh) is deliberately NOT checked. It targets
wasm32, where every .S here is guarded out to an empty object and fp.c takes
the portable path, so listing them would be noise rather than coverage.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CONSUMERS = [
    ("CMakeLists.txt",                  ROOT / "CMakeLists.txt"),
    ("bindings/python/build_ffi.py",    ROOT / "bindings" / "python" / "build_ffi.py"),
]


def main() -> int:
    asm = sorted(p.name for p in (ROOT / "src" / "arith").glob("*.S"))
    if not asm:
        print("error: no .S files under src/arith -- this check is looking in "
              "the wrong place", file=sys.stderr)
        return 2

    missing = []
    for label, path in CONSUMERS:
        if not path.is_file():
            print(f"error: {label} not found", file=sys.stderr)
            return 2
        text = path.read_text()
        for name in asm:
            if not re.search(re.escape(name), text):
                missing.append((label, name))

    print(f"{len(asm)} assembly sources, {len(CONSUMERS)} build files: "
          f"{len(missing)} missing")
    for label, name in missing:
        print(f"  {label} does not list {name}")

    if missing:
        print("\nA build that does not compile an .S file still links, right up "
              "until something calls into it. Add the file to the list above.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
