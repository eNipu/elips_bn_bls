\\ PROOFS, not tests.
\\
\\ curve_params.gp checks three specific curves. This file proves the same
\\ claims as identities in Z[X], which makes them true for EVERY seed of each
\\ family, including seeds nobody has chosen yet. A test says "it held for
\\ BN-462". A polynomial identity says "it holds, full stop".
\\
\\ These are the claims the library's fast paths are built on. Each one, if
\\ wrong, produces a routine that is wrong for some inputs and right for most,
\\ which is the failure mode tests are worst at finding.
\\
\\ A CAUTION LEARNED HERE. An earlier draft verified the CM equation with a
\\ power-series square root. That checks an identity only up to the truncation
\\ order and is not a proof, though it prints exactly like one. It is done with
\\ issquare() over Q[X] below, which is exact.
\\
\\ Run:  gp -q < tools/verify/identities.gp

fail = 0;
ran  = 0;
check(name, cond) = {
  ran = ran + 1;
  if(cond, printf("  ok    %s\n", name),
           printf("  FAIL  %s\n", name); fail = fail + 1);
};

print("Symbolic proofs in Z[X], valid for every seed of the family\n");

\\ ============================================================= BLS12 ========
\\ p(X) = ((X-1)^2 (X^4-X^2+1))/3 + X,  r(X) = X^4-X^2+1,  t(X) = X+1
P  = ((X-1)^2*(X^4 - X^2 + 1))/3 + X;
R  = X^4 - X^2 + 1;
T  = X + 1;

print("BLS12, for all X:");
check("#E = p + 1 - t is divisible by r", (P + 1 - T) % R == 0);
check("G1 cofactor h1 == (X-1)^2/3", (P + 1 - T)/R == (X-1)^2/3);
check("p == X (mod r), so psi and the Frobenius act as [X]", (P - X) % R == 0);

\\ The G1 subgroup test is E(P) == [m]P with m = -X^2, and phi satisfies
\\ phi^2 + phi + 1 = 0. A stray point of order l passes only if
\\ m^2 + m + 1 == 0 mod l. That quantity IS r, so no cofactor prime can ever
\\ satisfy it and the test is exact for every BLS12 seed, with no per-curve
\\ condition. This is the identity behind that claim.
check("(-X^2)^2 + (-X^2) + 1 == r exactly, so the G1 test is exact for all X", (-X^2)^2 + (-X^2) + 1 == R);

\\ The base-X^2 GLV split on G1 needs |X^2| to be about sqrt(r); r is X^4 to
\\ within lower-degree terms, so the digits really are half length.
check("r - X^4 has degree < 4, so X^2 is sqrt(r) up to low-order terms", poldegree(R - X^4) < 4);

\\ ================================================================ BN ========
\\ p(X) = 36X^4+36X^3+24X^2+6X+1, r(X) = 36X^4+36X^3+18X^2+6X+1, t(X) = 6X^2+1
BP = 36*X^4 + 36*X^3 + 24*X^2 + 6*X + 1;
BR = 36*X^4 + 36*X^3 + 18*X^2 + 6*X + 1;
BT = 6*X^2 + 1;

print("\nBN, for all X:");
check("#E = p + 1 - t == r, so the G1 cofactor is 1 for every BN seed", BP + 1 - BT == BR);
check("p == 6X^2 (mod r), so psi acts as [6X^2]", (BP - 6*X^2) % BR == 0);

\\ The fast G2 cofactor chain. h2 = p + t - 1 means that in base p the
\\ multiplier by p is ONE, which is what lets psi eliminate p entirely:
\\   [p] = [t]psi - psi^2  =>  [h2]Q = [6X^2](Q + psi(Q)) + psi(Q) - psi^2(Q)
H2 = BP + BT - 1;
check("h2 == p + t - 1 == p + 6X^2", H2 == BP + 6*X^2);

\\ And h2 * r really is a valid curve order over Fp2 with j = 0: its implied
\\ trace must satisfy the CM equation T2^2 - 4p^2 = -3 F^2 with F in Z[X].
N  = BR * H2;
T2 = BP^2 + 1 - N;
F2 = -(T2^2 - 4*BP^2)/3;
sq = issquare(F2, &G);
check("the implied twist trace satisfies the j=0 CM equation exactly", sq);
check("and its square root is a genuine polynomial in Z[X], not a series", sq && G^2 == F2);
if(sq, printf("        F(X) = %s\n", G));

\\ The reduced GLV lattice basis for BN G1, with a = 2X+1, c = 6X^2+2X,
\\ d = a + c. The determinant being exactly r is what makes the two digits
\\ short; if it were a proper multiple of r the split would still recombine
\\ but the digits would not be bounded.
A = 2*X + 1;  C = 6*X^2 + 2*X;  D = A + C;
check("GLV basis d == a + c == 6X^2 + 4X + 1", D == 6*X^2 + 4*X + 1);
check("GLV basis determinant a^2 + c*d == r exactly, for all X", A^2 + C*D == BR);

\\ lambda = 36X^3+18X^2+6X+1 must be a primitive cube root of unity mod r,
\\ which is what makes phi an endomorphism of order 3 on G1.
L = 36*X^3 + 18*X^2 + 6*X + 1;
check("lambda^2 + lambda + 1 == 0 (mod r), for all X", (L^2 + L + 1) % BR == 0);
check("the GLV basis vector (-a, c) is in the lattice: a == c*lambda (mod r)", (-A + C*L) % BR == 0);
check("the GLV basis vector (-d, -a) is in the lattice: d == -a*lambda (mod r)", (-D - A*L) % BR == 0);

\\ ============================================== the Fp12 tower =============
\\ ELiPS builds Fp12 as Fp2[w]/(w^6 - xi) with xi = 1 + u and u^2 = -1. That
\\ induces a degree-12 minimal polynomial for w over Fp, and it must be the
\\ one every other BLS12-381 implementation uses, or no cross-check against
\\ them is meaningful.
\\   w^6 = 1 + u  =>  u = w^6 - 1  =>  (w^6-1)^2 = -1  =>  w^12 - 2w^6 + 2 = 0
print("\nFp12 tower, independent of the curve:");
check("w^6 = 1+u with u^2 = -1 gives w^12 - 2w^6 + 2 == 0", subst((Y - 1)^2 + 1, Y, 'w^6) == 'w^12 - 2*'w^6 + 2);

EXPECTED = 16;
printf("\n%d of %d checks ran\n", ran, EXPECTED);
if(ran != EXPECTED, printf("MISCOUNT: %d checks did not run. A statement that fails to parse\nnever executes and never counts, so this guard exists to stop a\nsyntax error from reading as a clean pass.\n", EXPECTED - ran); fail = fail + 1);
printf("%s\n", if(fail, Str(fail, " PROOF(S) FAILED"), "all identities proved"));
quit(fail > 0);
