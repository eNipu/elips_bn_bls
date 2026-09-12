"""
Builds the CFFI extension that wraps include/elips/bls.h.

The C sources are compiled directly here rather than linked against a
pre-built library, so `pip install .` works from a checkout with nothing but a
C compiler: no CMake step, and no GMP, because as of issue #24 the library has
no runtime dependencies at all. That is the whole reason this package can be a
self-contained wheel.
"""
import os
import platform
import sys

from cffi import FFI


def _teach_distutils_about_assembly():
    """distutils only recognises .c, .cc, .cpp and friends as sources.

    Without this the .S files below are rejected outright, and if they were
    merely skipped instead the wheel would be 1.36x to 1.77x slower with
    nothing to say why. Silently slower is the worst outcome available, so
    this returns whether it worked and the caller acts on the answer.
    """
    try:
        from distutils.unixccompiler import UnixCCompiler
    except Exception:
        return False
    for ext in (".S", ".s"):
        if ext not in UnixCCompiler.src_extensions:
            UnixCCompiler.src_extensions.append(ext)
    return True


ASM_ENABLED = _teach_distutils_about_assembly()

HERE = os.path.dirname(os.path.abspath(__file__))
VENDOR = os.path.join(HERE, "vendor")
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

# Where the C comes from, and why it is copied rather than referenced.
#
# distutils refuses an absolute source path outright, and a relative one
# pointing at ../../src escapes the directory pip copies when it builds a
# source distribution, so the sdist would build only inside a git checkout.
# Copying the two directories the extension needs into vendor/ makes the
# package self-contained: the same tree builds from a checkout, from an sdist
# and into a wheel.
#
# The copy is refreshed on every build from a checkout, so it cannot go stale
# against the C it is meant to wrap.


def _sync_vendor():
    import shutil
    if not os.path.isdir(os.path.join(REPO, "src", "bls")):
        if os.path.isdir(os.path.join(VENDOR, "src", "bls")):
            return                                  # building from an sdist
        raise SystemExit("elips: no C sources found in %s or %s" % (REPO, VENDOR))
    for sub in ("include", "src"):
        dst = os.path.join(VENDOR, sub)
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(os.path.join(REPO, sub), dst)


_sync_vendor()


def _rel(*parts):
    """Relative to setup.py's directory: distutils rejects absolute paths."""
    return os.path.relpath(os.path.join(VENDOR, *parts), HERE)


C_SOURCES = [
    _rel("src", "arith", "fp.c"),
    _rel("src", "arith", "fpx.c"),
    _rel("src", "arith", "ec.c"),
    _rel("src", "arith", "serialize.c"),
    _rel("src", "arith", "wide.c"),
    _rel("src", "pairing", "miller.c"),
    _rel("src", "hash", "sha256.c"),
    _rel("src", "hash", "hash_to_curve.c"),
    _rel("src", "util", "random.c"),
    _rel("src", "util", "sysrand.c"),
    _rel("src", "bls", "bls.c"),
    _rel("src", "bls", "hkdf.c"),
]

# The assembly Montgomery multiply, worth 1.36x to 1.77x on a whole pairing
# (issue #7). Each file is guarded internally by #if on the target
# architecture, so listing both on any machine compiles one of them to an
# empty object rather than failing. distutils does not know .S is a source
# extension; setup.py teaches it, and drops these if that does not take.
ASM_SOURCES = [
    _rel("src", "arith", "fp_mul_x86_64.S"),
    _rel("src", "arith", "fp_mul_aarch64.S"),
]

# Exactly the public surface of include/elips/bls.h. Written out rather than
# fed the header, because the header carries an #error guard and macros cffi
# does not parse. If the two ever drift, the C compiler says so at build time:
# a mismatched prototype is an error, not a silent wrong call.
CDEF = """
int elips_bls_keygen(uint8_t *sk, const uint8_t *ikm, size_t ikm_len);
int elips_bls_keygen_random(uint8_t *sk);
int elips_bls_sk_to_pk(uint8_t *pk, const uint8_t *sk);
int elips_bls_pk_validate(const uint8_t *pk);

int elips_bls_sign(uint8_t *sig, const uint8_t *sk,
                   const uint8_t *msg, size_t msg_len);
int elips_bls_verify(const uint8_t *pk, const uint8_t *msg, size_t msg_len,
                     const uint8_t *sig);

int elips_bls_aggregate(uint8_t *out, const uint8_t *sigs, size_t n);
int elips_bls_aggregate_verify(const uint8_t *pks, size_t n,
                               const uint8_t *const *msgs, const size_t *msg_lens,
                               const uint8_t *sig);
int elips_bls_fast_aggregate_verify(const uint8_t *pks, const uint8_t *pops,
                                    size_t n, const uint8_t *msg, size_t msg_len,
                                    const uint8_t *sig);

int elips_bls_pop_prove(uint8_t *proof, const uint8_t *sk);
int elips_bls_pop_verify(const uint8_t *pk, const uint8_t *proof);

const char *elips_bls_strerror(int code);
"""


def make_ffi(with_asm=True):
    ffi = FFI()
    ffi.cdef(CDEF)
    sources = list(C_SOURCES) + (list(ASM_SOURCES) if with_asm else [])
    extra = ["-std=c11", "-O2"]
    if platform.system() != "Windows":
        extra.append("-fvisibility=hidden")
    ffi.set_source(
        "elips._elips",
        '#include "elips/bls.h"',
        sources=sources,
        include_dirs=[_rel("include"), _rel("src"), _rel("src", "bls")],
        define_macros=[("ELIPS_CURVE_BLS12_381", "1")],
        extra_compile_args=extra,
    )
    return ffi


if not ASM_ENABLED:
    print("elips: building without the assembly backend; the portable C will be "
          "used and this build will be slower. See issue #7.", file=sys.stderr)

ffibuilder = make_ffi(with_asm=ASM_ENABLED)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True, tmpdir=sys.argv[1] if len(sys.argv) > 1 else ".")
