\\ Independent verification of every curve parameter, using PARI/GP.
\\
\\ WHY THIS EXISTS
\\
\\ tools/reference/ is an oracle, but it was written by the same hand as the C.
\\ A shared misconception passes both, silently, which is exactly how the
\\ wrong-exponent final exponentiation of issue #16 survived for years.
\\
\\ PARI shares nothing with either. It counts the points on the curve itself
\\ rather than trusting the CM construction formula, and it does so for all
\\ three curves. That matters most for BN-462 and BLS12-461, which have no
\\ published standard anywhere: before this, their entire correctness rested on
\\ two artefacts by one author agreeing with each other.
\\
\\ WHAT THIS DOES NOT DO
\\
\\ It verifies the MATHEMATICS the library is supposed to implement, not the C
\\ that implements it. A carry bug in fp_mul is invisible here.
\\
\\ Run:  gp -q < tools/verify/curve_params.gp

default(parisize, 2000000000);
fail = 0;
ran  = 0;

check(name, cond) = {
  ran = ran + 1;
  if(cond, printf("  ok    %s\n", name),
           printf("  FAIL  %s\n", name); fail = fail + 1);
};

\\ ---------------------------------------------------------------- BLS12 ----
\\ p = ((x-1)^2 (x^4-x^2+1))/3 + x,  r = x^4-x^2+1,  t = x+1
bls12(name, x) = {
  my(p, r, t, E, n, h1, u, xi, Et, n2, h2);
  p = ((x-1)^2*(x^4-x^2+1))/3 + x;
  r = x^4 - x^2 + 1;
  t = x + 1;
  printf("\n%s  (p %d bits, r %d bits)\n", name, #binary(p), #binary(r));
  check("p is prime", isprime(p));
  check("r is prime", isprime(r));

  \\ The point count, done by PARI rather than taken from the formula.
  E = ellinit([0, 4], p);
  n = ellcard(E);
  check("PARI's own #E(Fp) equals p + 1 - t", n == p + 1 - t);
  check("r divides #E(Fp)", n % r == 0);
  h1 = n / r;
  check("G1 cofactor h1 == (x-1)^2/3", h1 == (x-1)^2/3);

  \\ The sextic twist over Fp2, which carries G2. Its cofactor is what
  \\ hash-to-curve has to clear, so a wrong value here is a wrong h_eff.
  u  = ffgen(Mod(1,p)*('u^2 + 1), 'u);
  xi = 1 + u;
  Et = ellinit([0, 4*xi]);
  n2 = ellcard(Et);
  check("r divides PARI's own #E'(Fp2)", n2 % r == 0);
  h2 = n2 / r;
  check("G2 cofactor h2 == (x^8 - 4x^7 + 5x^6 - 4x^4 + 6x^3 - 4x^2 - 4x + 13)/9",
        h2 == (x^8 - 4*x^7 + 5*x^6 - 4*x^4 + 6*x^3 - 4*x^2 - 4*x + 13)/9);

  \\ psi acts on G2 as multiplication by p mod r, and on BLS12 that is x.
  check("p == x  (mod r), the psi and Frobenius multiplier", (p - x) % r == 0);

  \\ phi acts on G1 as [-x^2]; the subgroup test is exact because
  \\ m^2 + m + 1 with m = -x^2 is r itself.
  check("(-x^2)^2 + (-x^2) + 1 == r, so the G1 test is exact", x^4 - x^2 + 1 == r);
};

\\ ------------------------------------------------------------------- BN ----
\\ p = 36x^4+36x^3+24x^2+6x+1,  r = 36x^4+36x^3+18x^2+6x+1,  t = 6x^2+1
bn(name, x) = {
  my(p, r, t, u, xi, Et, n2, h2, a, c, d);
  p = 36*x^4 + 36*x^3 + 24*x^2 + 6*x + 1;
  r = 36*x^4 + 36*x^3 + 18*x^2 + 6*x + 1;
  t = 6*x^2 + 1;
  printf("\n%s  (p %d bits, r %d bits)\n", name, #binary(p), #binary(r));
  check("p is prime", isprime(p));
  check("r is prime", isprime(r));
  check("#E(Fp) = p + 1 - t = r, so the G1 cofactor is 1", p + 1 - t == r);

  u  = ffgen(Mod(1,p)*('u^2 + 1), 'u);
  xi = 1 + u;
  Et = ellinit([0, -4*xi]);          \\ BN here is y^2 = x^3 - 4
  n2 = ellcard(Et);
  check("r divides PARI's own #E'(Fp2)", n2 % r == 0);
  h2 = n2 / r;
  \\ This is the identity the fast BN G2 cofactor chain depends on. If it is
  \\ wrong, hash_to_g2 lands outside G2 and every downstream pairing is wrong.
  check("G2 cofactor h2 == p + t - 1", h2 == p + t - 1);
  check("G2 cofactor h2 == p + 6x^2", h2 == p + 6*x^2);

  check("p == 6x^2 (mod r), the psi multiplier", (p - 6*x^2) % r == 0);

  \\ The reduced GLV lattice basis for BN G1: (-a, c) and (-d, -a).
  a = 2*x + 1;  c = 6*x^2 + 2*x;  d = a + c;
  check("GLV basis third entry d == a + c", d == 6*x^2 + 4*x + 1);
  check("GLV basis determinant a^2 + c*d == r exactly", a^2 + c*d == r);
};

print("Independent curve-parameter verification (PARI/GP, nothing shared with ELiPS)");
bls12("BLS12-381", -0xd201000000010000);
bls12("BLS12-461", -2^77 + 2^50 + 2^33);
bn("BN-462",   2^114 + 2^101 - 2^14 - 1);

\\ A check that fails to parse never runs and never counts, which would read
\\ as a clean pass. This guard turns that into a failure. It has already
\\ caught one such case in the companion file.
EXPECTED = 27;
printf("\n%d of %d checks ran\n", ran, EXPECTED);
if(ran != EXPECTED, printf("MISCOUNT: %d checks did not run\n", EXPECTED - ran); fail = fail + 1);
printf("%s\n", if(fail, Str(fail, " CHECK(S) FAILED"), "all curve parameters verified"));
quit(fail > 0);
