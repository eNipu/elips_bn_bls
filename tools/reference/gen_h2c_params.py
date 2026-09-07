"""Emit the hash-to-curve constants for each curve (Phase 5b, issue #6).

Kept separate from gen_params.py because the isogeny tables are large and only
one curve has them; mixing them into the field parameters would bury the part a
reader actually needs to check.

Nothing is emitted before tools/reference/h2c_ref.py's self-test passes, which
is what proves the isogeny data really is an isogeny and that the whole chain
reproduces the RFC 9380 test vectors. A wrong constant here produces points that
are on the curve and in the group and simply disagree with every other
implementation, so "it runs" is not evidence.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from elips_ref import CURVES                                    # noqa: E402
import h2c_ref                                                  # noqa: E402

W = 64
MASK = (1 << W) - 1
SEP = ",\n        "


def limbs(x, n):
    return SEP.join("UINT64_C(0x%016x)" % ((x >> (W * i)) & MASK) for i in range(n))


def emit_curve(name, macro):
    cv = CURVES[name]
    p = cv.p
    n = (p.bit_length() + W - 1) // W
    R = 1 << (W * n)

    def mont(v):
        return (v * R) % p

    def fp_const(nm, f):
        return ("  static const limb_t %s[%d] = {\n        %s\n  };\n"
                % (nm, n, limbs(mont(f.v), n)))

    def fp2_const(nm, f):
        a, b = f.limbs()
        return ("  static const limb_t %s[2][%d] = {\n    { %s },\n    { %s }\n  };\n"
                % (nm, n, limbs(mont(a), n), limbs(mont(b), n)))

    def fp_table(nm, polys):
        out = ""
        for tag, poly in zip(("XNUM", "XDEN", "YNUM", "YDEN"), polys):
            out += "  #define %s_%s_LEN %d\n" % (nm, tag, len(poly))
            rows = ",\n".join("    { %s }" % limbs(mont(c.v), n) for c in poly)
            out += ("  static const limb_t %s_%s[%d][%d] = {\n%s\n  };\n"
                    % (nm, tag, len(poly), n, rows))
        return out

    def fp2_table(nm, polys):
        out = ""
        for tag, poly in zip(("XNUM", "XDEN", "YNUM", "YDEN"), polys):
            out += "  #define %s_%s_LEN %d\n" % (nm, tag, len(poly))
            rows = []
            for c in poly:
                a, b = c.limbs()
                rows.append("    { { %s },\n      { %s } }"
                            % (limbs(mont(a), n), limbs(mont(b), n)))
            out += ("  static const limb_t %s_%s[%d][2][%d] = {\n%s\n  };\n"
                    % (nm, tag, len(poly), n, ",\n".join(rows)))
        return out

    g1 = h2c_ref.suite(cv, "G1")
    g2 = h2c_ref.suite(cv, "G2")
    assert g1.kind == g2.kind

    txt = "#if defined(%s)\n" % macro
    txt += '  #define ELIPS_H2C_SUITE_G1 "%s"\n' % g1.suite_id
    txt += '  #define ELIPS_H2C_SUITE_G2 "%s"\n' % g2.suite_id
    txt += "  /* RFC 9380 5.3: ceil((ceil(log2(p)) + 128) / 8). */\n"
    txt += "  #define ELIPS_H2C_L        %d\n" % g1.L

    if g1.kind == "SSWU":
        txt += ("  #define ELIPS_H2C_SSWU     1\n"
                "  /* The map runs on a curve isogenous to this one, because\n"
                "   * simplified SWU needs A*B != 0 and these curves have A = 0. */\n")
        txt += fp_const("H2C_G1_ISO_A", g1.iso_A)
        txt += fp_const("H2C_G1_ISO_B", g1.iso_B)
        txt += fp_const("H2C_G1_Z", g1.Z)
        # -B/A and B/(Z*A) are the only two quotients the map needs, and both
        # are functions of the constants above. Precomputing them keeps two
        # constant-time inversions -- 26 us each -- out of every call.
        txt += fp_const("H2C_G1_MB_OVER_A", -g1.iso_B * g1.iso_A.inv0())
        txt += fp_const("H2C_G1_B_OVER_ZA", g1.iso_B * (g1.Z * g1.iso_A).inv0())
        txt += fp_table("H2C_G1_ISO", g1.iso)
        txt += fp2_const("H2C_G2_ISO_A", g2.iso_A)
        txt += fp2_const("H2C_G2_ISO_B", g2.iso_B)
        txt += fp2_const("H2C_G2_Z", g2.Z)
        txt += fp2_const("H2C_G2_MB_OVER_A", -g2.iso_B * g2.iso_A.inv0())
        txt += fp2_const("H2C_G2_B_OVER_ZA", g2.iso_B * (g2.Z * g2.iso_A).inv0())
        txt += fp2_table("H2C_G2_ISO", g2.iso)
    else:
        txt += ("  #define ELIPS_H2C_SVDW     1\n"
                "  /* Shallue-van de Woestijne needs no isogeny; c1..c4 are the\n"
                "   * precomputed values of RFC 9380 6.6.1. Z was chosen by the\n"
                "   * Appendix H.1 search, not picked by hand. */\n")
        txt += "  #define ELIPS_H2C_G1_Z_INT %d\n" % g1.z_int
        txt += fp_const("H2C_G1_Z", g1.Z)
        for i, c in enumerate(g1.C):
            txt += fp_const("H2C_G1_C%d" % (i + 1), c)
        txt += "  #define ELIPS_H2C_G2_Z_INT %d\n" % g2.z_int
        txt += fp2_const("H2C_G2_Z", g2.Z)
        for i, c in enumerate(g2.C):
            txt += fp2_const("H2C_G2_C%d" % (i + 1), c)

    # The fast G2 cofactor chain, where one exists and has been verified.
    fast = h2c_ref.g2_fast_clear_coeffs(cv)
    if fast:
        a, b = fast
        txt += ("  /* Budroni-Pintore (ePrint 2017/419):\n"
                "   *   [h_eff]Q = [x^2-x-1]Q + [x-1]psi(Q) + psi^2([2]Q)\n"
                "   * Two ladders over the short parameter replace one over h_eff.\n"
                "   * Verified against [h_eff] on random points of the twist before\n"
                "   * these constants were emitted. */\n"
                "  #define ELIPS_H2C_G2_FAST_CLEAR 1\n")
        for nm, v in (("A", a), ("B", b)):
            av = abs(v)
            an = (av.bit_length() + W - 1) // W
            txt += "  #define H2C_G2_CLEAR_%s_BITS %d\n" % (nm, av.bit_length())
            txt += "  #define H2C_G2_CLEAR_%s_NEG  %d\n" % (nm, 1 if v < 0 else 0)
            txt += ("  static const limb_t H2C_G2_CLEAR_%s[%d] = {\n        %s\n  };\n"
                    % (nm, an, limbs(av, an)))

    txt += "  /* clear_cofactor multipliers. */\n"
    for tag, suite in (("G1", g1), ("G2", g2)):
        h = suite.h_eff
        hn = (h.bit_length() + W - 1) // W
        txt += "  #define H2C_HEFF_%s_BITS  %d\n" % (tag, h.bit_length())
        txt += ("  static const limb_t H2C_HEFF_%s[%d] = {\n        %s\n  };\n"
                % (tag, hn, limbs(h, hn)))

    txt += "#endif\n"
    return txt


HEADER = """/* Generated by tools/reference/gen_h2c_params.py -- do not edit.
 *
 * RFC 9380 hash-to-curve constants. Field elements are in Montgomery form, so
 * the library needs no initialisation step for them.
 *
 * BLS12-381 uses simplified SWU over an isogenous curve, which is the suite
 * RFC 9380 registers, so its output matches every conforming implementation.
 * BLS12-461 and BN-462 have no registered suite and use Shallue-van de
 * Woestijne, which RFC 9380 defines for exactly this case (A = 0, no isogeny
 * available) and uses itself for BN254.
 *
 * Everything here except the BLS12-381 isogeny tables is derived from the curve
 * parameters. The isogeny tables come from tools/reference/h2c_iso_bls12_381.json
 * and are verified to be an isogeny -- and to reproduce the RFC's own test
 * vectors -- before this file is written.
 */
#ifndef ELIPS_H2C_PARAMS_H
#define ELIPS_H2C_PARAMS_H

#include "elips/fp_params.h"

"""

FOOTER = """
#ifndef ELIPS_H2C_L
#  error "No curve selected, or this curve has no hash-to-curve parameters."
#endif

#endif /* ELIPS_H2C_PARAMS_H */
"""

if __name__ == "__main__":
    if h2c_ref.selftest(verbose=False):
        sys.exit("h2c_ref self-test failed; refusing to emit constants")

    parts = [HEADER]
    for name, macro in (("BLS12-381", "ELIPS_CURVE_BLS12_381"),
                        ("BLS12-461", "ELIPS_CURVE_BLS12_461"),
                        ("BN-462",    "ELIPS_CURVE_BN_462")):
        parts.append(emit_curve(name, macro))
    parts.append(FOOTER)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "include", "elips", "h2c_params.h")
    open(out, "w").write("\n".join(parts))
    print("wrote include/elips/h2c_params.h")
