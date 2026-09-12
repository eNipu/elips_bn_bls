#!/bin/sh
# Builds the WebAssembly module from the C in ../../src.
#
#   . /path/to/emsdk/emsdk_env.sh && ./build.sh
#
# No CMake and no GMP: as of issue #24 the library has no runtime dependencies,
# which is the only reason a browser build is a build script rather than a
# project. The .S files are listed and compile to nothing off their target, so
# the same command works on any host.
set -eu

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT="$HERE/dist"
mkdir -p "$OUT"

SRC="
  $ROOT/src/arith/fp.c
  $ROOT/src/arith/fpx.c
  $ROOT/src/arith/ec.c
  $ROOT/src/arith/serialize.c
  $ROOT/src/arith/wide.c
  $ROOT/src/pairing/miller.c
  $ROOT/src/hash/sha256.c
  $ROOT/src/hash/hash_to_curve.c
  $ROOT/src/util/random.c
  $ROOT/src/util/sysrand.c
  $ROOT/src/bls/bls.c
  $ROOT/src/bls/hkdf.c
"

# Only the BLS surface is exported. Nothing else needs to cross into
# JavaScript, and a smaller export list is a smaller module.
FUNCS='["_elips_bls_keygen","_elips_bls_keygen_random","_elips_bls_sk_to_pk",
        "_elips_bls_pk_validate","_elips_bls_sign","_elips_bls_verify",
        "_elips_bls_aggregate","_elips_bls_aggregate_verify",
        "_elips_bls_fast_aggregate_verify","_elips_bls_pop_prove",
        "_elips_bls_pop_verify","_elips_bls_strerror","_malloc","_free"]'

emcc $SRC \
  -O3 -std=c11 \
  -DELIPS_CURVE_BLS12_381 \
  -I"$ROOT/include" -I"$ROOT/src" -I"$ROOT/src/bls" \
  -s EXPORTED_FUNCTIONS="$FUNCS" \
  -s EXPORTED_RUNTIME_METHODS='["HEAPU8","UTF8ToString"]' \
  -s MODULARIZE=1 \
  -s EXPORT_ES6=1 \
  -s ENVIRONMENT='web,worker,node' \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s SINGLE_FILE=1 \
  -s FILESYSTEM=0 \
  -o "$OUT/elips_wasm.mjs"

echo "built $OUT/elips_wasm.mjs  ($(du -h "$OUT/elips_wasm.mjs" | cut -f1))"
