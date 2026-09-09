# Understanding and optimizing a BLS12-381 Miller loop

## A guided research note, from field arithmetic to working code

**Audience:** an entry-level researcher who knows basic programming and algebra.
**Written:** 2026-09-10. **Measurements discussed:** 2026-09-09.
**Scope:** the preserved ELiPS research prototypes, not a production patch.

The main idea is simple:

> An intermediate result does not always need to be exact. It needs to
> preserve everything that the final computation observes.

In a pairing, final exponentiation removes certain multiplicative factors.
We can sometimes keep those factors instead of spending operations to remove
them earlier.

This note develops two examples:

1. A prepared-line multiplication using **8 rather than 15 Fp2
   multiplications** in the compared implementations.
2. A general squaring using **30 rather than 36 base-field
   multiplications**, by returning four times the square.

Both preserve the reduced pairing, but **neither preserves every raw
intermediate value**. That distinction is the central correctness condition.

These are tested implementation improvements assembled from established
ideas. **Publication-level novelty and optimality are not established.**
Earlier overly strong claims are corrected in
[CORRECTIONS.md](CORRECTIONS.md); the original experiments remain preserved.

### Contents

1. [Learning path](#1-learning-path)
2. [Notation and cost model](#2-notation-and-cost-model)
3. [Finite fields and the tower](#3-finite-fields-and-the-tower)
4. [Pairings and the Miller loop](#4-pairings-and-the-miller-loop)
5. [Why a scale factor can disappear](#5-why-a-scale-factor-can-disappear)
6. [Normalize a prepared line](#6-normalize-a-prepared-line)
7. [Derive the 8M2 line kernel](#7-derive-the-8m2-line-kernel)
8. [Derive Fourier-based squaring](#8-derive-fourier-based-squaring)
9. [Two hand-worked examples](#9-two-hand-worked-examples)
10. [Combine the kernels safely](#10-combine-the-kernels-safely)
11. [Read the code in order](#11-read-the-code-in-order)
12. [Run and reproduce the experiments](#12-run-and-reproduce-the-experiments)
13. [Interpret the measurements](#13-interpret-the-measurements)
14. [Correctness and security checklist](#14-correctness-and-security-checklist)
15. [Other research paths](#15-other-research-paths)
16. [Exercises and next experiments](#16-exercises-and-next-experiments)
17. [Further reading](#17-further-reading)

## 1. Learning path

You do not need to understand all of pairing theory before checking these
arithmetic identities. There are three separate questions:

1. **Local algebra:** does a kernel return exactly the formula it promises?
2. **Pairing correctness:** does the deliberate scale factor disappear at
   the end?
3. **Performance:** does the resulting program actually run faster?

Prove or test them separately. A good benchmark does not prove correctness,
and a correct formula is not automatically fast.

A suggested first pass:

- Read sections 2–5 for the mathematical setting.
- Work through either example in section 9 by hand.
- Read the corresponding derivation in section 7 or 8.
- Run [tutorial_examples.py](tutorial_examples.py).
- Read sections 13–14 before interpreting a timing or changing an API.

The companion program uses only Python's standard library and the existing
reference code. Its examples use public, deterministic values. It is not a
cryptographic application.

## 2. Notation and cost model

| Symbol | Meaning |
| --- | --- |
| \(p\) | Prime defining the base field; 381 bits for BLS12-381 |
| \(r\) | Prime order of the pairing groups; 255 bits here |
| \(\mathbb F_{p^d}\) | Field with \(p^d\) elements |
| \(\mathbb F_{p^d}^{\times}\) | Its nonzero elements under multiplication |
| \(P\in G_1,\ Q\in G_2\) | Input elliptic-curve points |
| \([m]Q\) | Scalar multiplication of a point, not field multiplication |
| \(f\) | Miller accumulator, an element of \(\mathbb F_{p^{12}}\) |
| \(L\) | A line function evaluated at the pairing's other argument |
| \(E=(p^{12}-1)/r\) | Final exponent |
| \(M_p,\ S_p,\ I_p\) | Base-field multiplication, squaring, inversion |
| \(M_d,\ S_d,\ I_d\) | Corresponding operation in \(\mathbb F_{p^d}\) |
| \(A_d\) | An addition or subtraction in \(\mathbb F_{p^d}\) |

An expression such as \(8M_2\) counts eight *general multiplications* in
\(\mathbb F_{p^2}\). It is not a time measurement and does not include all
additions, memory traffic, or setup work.

In the C routines used here:

$$
M_2=3M_p,\qquad S_2=2M_p,\qquad
M_6=6M_2.
$$

These equalities describe multiplication counts in these particular
implementations. They are not universal lower bounds.

The Python oracle sometimes uses simpler schoolbook formulas instead of
the C algorithm. Its value is independence and readability. **Do not use
Python execution time, or a count of Python `*` operators, as the C cost model.**

## 3. Finite fields and the tower

### 3.1 Start with modular arithmetic

The field \(\mathbb F_p\) consists of integers modulo a prime \(p\).
Addition and multiplication wrap around modulo \(p\).

For example, in \(\mathbb F_{13}\):

$$
9+8=4,\qquad 7\cdot 8=4,\qquad 5^{-1}=8,
$$

because \(17\equiv4\), \(56\equiv4\), and \(5\cdot8\equiv1\pmod{13}\).
Division means multiplication by an inverse. **Zero has no inverse.**

The final example in section 9 uses a small field for hand calculation.
That small field is not a secure pairing curve.

### 3.2 Why extension fields appear

The pairing output is not a single base-field element. For BLS12-381 it
lives in a subgroup of \(\mathbb F_{p^{12}}^\times\).

An extension field is implemented with polynomial coefficients and a
reduction rule. It is not implemented by allocating a table containing
all \(p^{12}\) elements.

The repository uses:

$$
\begin{aligned}
\mathbb F_{p^2}&=\mathbb F_p[u]/(u^2+1),\\
\xi&=1+u,\\
\mathbb F_{p^6}&=\mathbb F_{p^2}[v]/(v^3-\xi),\\
\mathbb F_{p^{12}}&=\mathbb F_{p^6}[w]/(w^2-v).
\end{aligned}
$$

The notation means that calculations obey

$$
u^2=-1,\qquad v^3=\xi,\qquad w^2=v,\qquad w^6=\xi.
$$

The defining polynomials must be irreducible for these quotients to be
fields. That condition holds for the configured tower; it is not true for
an arbitrary choice of constants.

### 3.3 Fp2 arithmetic

An Fp2 element is \(a_0+a_1u\). Multiplication expands to

$$
(a_0+a_1u)(b_0+b_1u)
=(a_0b_0-a_1b_1)+(a_0b_1+a_1b_0)u.
$$

Computing all four products is unnecessary. Let

$$
t_0=a_0b_0,\quad t_1=a_1b_1,\quad
t_2=(a_0+a_1)(b_0+b_1).
$$

Then the result is

$$
(t_0-t_1)+(t_2-t_0-t_1)u.
$$

This is Karatsuba's idea: three multiplications plus extra additions.

Squaring is cheaper:

$$
(a_0+a_1u)^2=(a_0+a_1)(a_0-a_1)+2a_0a_1u.
$$

The C routine therefore uses two base-field multiplications.

Some constant multiplications are especially cheap:

$$
\begin{aligned}
u(a_0+a_1u)&=-a_1+a_0u,\\
\xi(a_0+a_1u)&=(a_0-a_1)+(a_0+a_1)u.
\end{aligned}
$$

These need swaps, negation, and additions, not a general Fp2 multiplication.
“Cheap” still does not mean “zero execution time.”

### 3.4 An Fp12 element has several useful views

The normal tower view is

$$
f=A+Dw,\qquad
A,D\in\mathbb F_{p^6}.
$$

Flatten it as

$$
f=g_0+g_1w+g_2w^2+g_3w^3+g_4w^4+g_5w^5,
\qquad g_i\in\mathbb F_{p^2}.
$$

The storage order is:

```text
C:      f[0] = (g0, g2, g4),     f[1] = (g1, g3, g5)
Python: f.d0 = Fp6(g0, g2, g4),  f.d1 = Fp6(g1, g3, g5)
```

For the squaring optimization, define \(s=w^3\). Then \(s^2=\xi\), so:

$$
\mathbb F_{p^4}=\mathbb F_{p^2}[s]/(s^2-\xi),\qquad
\mathbb F_{p^{12}}=\mathbb F_{p^4}[w]/(w^3-s).
$$

Now write

$$
f=a+bw+cw^2
$$

with

$$
a=g_0+g_3s,\quad b=g_1+g_4s,\quad c=g_2+g_5s.
$$

This is the **same field element**, regrouped into three Fp4 coordinates.
It is not a conversion into a new cryptosystem.

| Coordinate | Fp2 components | Python source |
| --- | --- | --- |
| \(a\) | \((g_0,g_3)\) | `(f.d0.c0, f.d1.c1)` |
| \(b\) | \((g_1,g_4)\) | `(f.d1.c0, f.d0.c2)` |
| \(c\) | \((g_2,g_5)\) | `(f.d0.c1, f.d1.c2)` |

The regrouping requires no field multiplication, although copying
coordinates can still have a machine cost.

### 3.5 A note on Montgomery representation

The C base-field values are stored in Montgomery form, effectively as
\(aR\bmod p\) for a fixed power-of-two radix \(R\). This makes repeated
modular multiplication efficient.

All equations in this note describe ordinary field values. The C
`fp_*` routines maintain the internal representation. Do not put the
integer `1` directly into a C field array and expect it to represent the
field identity; use `fp_set_one`. The Python oracle uses ordinary residues.

## 4. Pairings and the Miller loop

### 4.1 What a pairing computes

A pairing is a map

$$
e:G_1\times G_2\longrightarrow G_T
$$

with bilinearity:

$$
e([a]P,[b]Q)=e(P,Q)^{ab}.
$$

The input groups use elliptic-curve point addition. The target group uses
field multiplication. Here all three groups have prime order \(r\).
Nondegeneracy rules out the useless map that always returns one.

For BLS12-381, the base curve is \(y^2=x^3+4\). G2 is represented on a
sextic twist over Fp2, letting its point arithmetic use Fp2 coordinates
instead of arbitrary Fp12 coordinates.

The family parameter is

$$
x=-\mathtt{0xd201000000010000},
$$

and

$$
r=x^4-x^2+1,\qquad
p=\frac{(x-1)^2r}{3}+x.
$$

The *embedding degree* is the smallest positive \(k\) with
\(r\mid p^k-1\); here \(k=12\).

### 4.2 The two stages

Conceptually:

```text
P, Q
  |
  v
Miller loop: point updates + field accumulator updates
  |
  v
raw f in Fp12
  |
  v
final exponentiation: f^E, E = (p^12-1)/r
  |
  v
reduced pairing in GT
```

The optimal-ate construction determines the short scalar loop and any
correction lines. This note does not reprove that pairing construction.
It changes the arithmetic used to execute an already validated schedule.

### 4.3 Why lines occur

A line through two curve points, or a tangent at a point, can be written as

$$
\ell(X,Y)=(Y-y_T)-\lambda(X-x_T),
$$

where \(\lambda\) is the chord or tangent slope. Miller's algorithm combines
such rational functions to obtain a function with a prescribed divisor.

For a basic doubling step, the complete rational-function update is
schematically

$$
f\leftarrow f^2
\frac{\ell_{T,T}(P)}{v_{[2]T}(P)},\qquad T\leftarrow[2]T,
$$

where \(v_{[2]T}\) is the vertical line through \([2]T\). The fraction is
multiplied into \(f^2\); it is not added.

In the denominator-eliminated pairing path used here, the omitted vertical
factors lie in a subfield whose nonzero elements final exponentiation
annihilates. The field update then looks like

$$
f\leftarrow f^2L_{\rm dbl},
$$

followed, on nonzero digits, by

$$
f\leftarrow fL_{\rm add}.
$$

This schematic form omits signed-digit initialization and other
construction details. Copy the tested schedule, not just this pseudocode,
when implementing a pairing.

### 4.4 Fixed-argument preparation and multi-pairing

If Q is reused, prepare all Q-dependent line coefficients once. Online
work then evaluates these lines at a new P and updates the accumulator.
This is a **prepared pairing**, not an unprepared one.

For several pairs, a shared iteration is

$$
f\leftarrow f^2\prod_j L_{j,\rm dbl}
$$

with any required addition lines also multiplied in. The accumulator is
squared once per iteration, not once per pair.

The current BLS12-381 signed-digit schedule uses **64 squarings and 69
lines per term**. These are properties of this schedule, not lower bounds.
Its absolute binary seed has 64 bits and permits a 63-doubling schedule;
changing schedules requires its own correctness and cost comparison.

## 5. Why a scale factor can disappear

### 5.1 The basic lemma

Let \(d\) divide 12, and suppose

$$
p^d-1\mid E.
$$

Every nonzero \(\alpha\in\mathbb F_{p^d}\) satisfies

$$
\alpha^{p^d-1}=1.
$$

Writing \(E=k(p^d-1)\), we have

$$
\alpha^E=(\alpha^{p^d-1})^k=1.
$$

Consequently, for nonzero \(f\),

$$
(\alpha f)^E=\alpha^Ef^E=f^E.
$$

**The raw values differ, but the final pairing values agree.**

The divisibility condition is essential. Merely calling something a
“subfield factor” is not a proof for an arbitrary curve or exponent.
The scripts check it for \(d=1,2,4,6\) on the three configured curves.

For example, because \(r\mid p^6+1\),

$$
E=(p^6-1)\frac{p^6+1}{r},
$$

so nonzero Fp6 factors disappear. For the line kernel we only need the
Fp2 case; for the new squaring we only need the Fp case.

### 5.2 View it as permitted information loss

It is useful to write

$$
f\sim g\quad\text{when}\quad
g=\alpha f,\ \alpha\in\mathbb F_{p^2}^{\times}.
$$

This is a sufficient relation for equal final pairings. It does not
describe every pair of field elements with equal final exponentiation.

The optimization searches for a cheaper representative of the same
permitted class. It does **not** approximate a field value numerically.
All arithmetic is exact.

### 5.3 Why repeated scaling remains safe

Suppose the optimized accumulator is \(\widetilde f=\alpha f\).
After squaring, the discrepancy is \(\alpha^2\).
Replacing a line \(L\) by \(\beta L\) gives discrepancy
\(\alpha^2\beta\). A subfield is closed under these operations, so the
discrepancy stays in the permitted nonzero subfield.

This induction is more useful than checking only one isolated kernel.
It explains why the transformations compose across a complete loop.

### 5.4 Where this reasoning stops

- It does not allow a zero scale.
- It does not preserve serialized raw Miller values.
- It does not authorize replacing the final pairing by an arbitrary value
  related by exponentiation.
- It does not preserve membership in every subgroup used by a specialized
  formula.

In particular, Granger–Scott squaring requires a cyclotomic input. The easy
final exponent \((p^6-1)(p^2+1)\) puts a nonzero value in that subgroup.
Raising only to \(p^6-1\) is not sufficient. Arbitrarily scaling a cyclotomic
value can also take it outside that subgroup.

## 6. Normalize a prepared line

### 6.1 The line shape

In this tower, a prepared line has the form

$$
L=\ell_a y_P+\ell_c w^3+\ell_b x_Pw^5,
\qquad \ell_a,\ell_b,\ell_c\in\mathbb F_{p^2}.
$$

The code calls these coefficients `a`, `b`, and `c`; the \(\ell\) prefix
here avoids confusion with the later Fp4 coordinates.

Only three of the six Fp2 coefficients are nonzero. This is a *sparse*
line. The production sparse kernel already uses this structure and costs
15M2 in the implementation under comparison.

### 6.2 Divide once in preparation, not once per line online

For a fixed Q, prepare

$$
\beta=\frac{\ell_c}{\ell_a},\qquad
\gamma=\frac{\ell_b}{\ell_a}.
$$

For each new P, compute

$$
\operatorname{iy}=y_P^{-1},\qquad
\operatorname{xy}=x_Py_P^{-1}.
$$

Then set

$$
B=\beta\operatorname{iy},\qquad
C=\gamma\operatorname{xy}.
$$

The normalized line is

$$
\widehat L=1+Bw^3+Cw^5,\qquad
L=(\ell_a y_P)\widehat L.
$$

We removed a nonzero Fp2 factor, which section 5 permits.

The leading coefficient being **one** is important: multiplying an
accumulator by that coefficient no longer costs a general multiplication.

### 6.3 Preconditions and costs

For a nonidentity point of odd prime order, \(y_P\ne0\), since a point with
zero y-coordinate has order two. The actual signed scalar prefixes are
also checked to ensure that each normalized line denominator is nonzero.
This does not justify normalizing arbitrary malformed or unvalidated data.

Batch inversion computes n inverses with one inversion and
\(3(n-1)\) multiplications:

1. Build prefix products \(a_0,\ a_0a_1,\ldots,a_0\cdots a_{n-1}\).
2. Invert the final product once.
3. Walk backward, extracting each inverse and updating the running inverse.

For 69 lines, normalizing both remaining coefficients costs

$$
1I_2+\bigl(3(69-1)+2\cdot69\bigr)M_2
=1I_2+342M_2.
$$

This is additional **offline** work. Online, each P costs an extra
\(1I_p+1M_p\). Each line still needs two Fp2-by-Fp scalings, costing \(4M_p\).

The normalized table stores two Fp2 values per line rather than three.
Its measured C size is 13,248 rather than 19,872 bytes.
Those sizes depend on this field representation; the one-third reduction
in coefficient count is the more general observation.

Preparation is worthwhile primarily when Q is reused. Do not hide its
cost when comparing one-use pairings.

## 7. Derive the 8M2 line kernel

### 7.1 Reduce the problem to two smaller polynomial products

Write \(f=A+Dw\), with \(A,D\in\mathbb F_{p^6}\). Since \(v=w^2\),

$$
\widehat L=1+w^3(B+Cv).
$$

Expanding gives

$$
f\widehat L
=\left[A+v^2D(B+Cv)\right]
 +\left[D+vA(B+Cv)\right]w.
$$

We need to multiply each of A and D by the same linear polynomial
\(B+Cv\).

### 7.2 Evaluate instead of computing every cross product

Let \(T(z)=t_0+t_1z+t_2z^2\) and consider

$$
R(z)=T(z)(B+Cz).
$$

Before reduction, write

$$
R(z)=r_0+r_1z+r_2z^2+r_3z^3.
$$

Four evaluations determine its four coefficients:

$$
\begin{aligned}
m_0&=t_0B,\\
m_3&=t_2C,\\
m_+&=(t_0+t_1+t_2)(B+C),\\
m_-&=(t_0-t_1+t_2)(B-C).
\end{aligned}
$$

Here “evaluation at infinity” means extracting the leading coefficient
\(r_3=t_2C\). It is not a numerical operation involving an infinite value.

Observe:

$$
\begin{aligned}
m_+&=r_0+r_1+r_2+r_3,\\
m_-&=r_0-r_1+r_2-r_3.
\end{aligned}
$$

Therefore

$$
2r_1=m_+-m_--2m_3,\qquad
2r_2=m_++m_--2m_0.
$$

Ordinary interpolation would divide these expressions by two.
Instead, keep the common factor two in **all** coefficients.

Reducing \(z^3=\xi\), define

$$
\begin{aligned}
H_0&=2(m_0+\xi m_3),\\
H_1&=m_+-m_--2m_3,\\
H_2&=m_++m_--2m_0.
\end{aligned}
$$

Then

$$
H(T)=H_0+H_1v+H_2v^2=2T(B+Cv).
$$

This costs four M2, not six schoolbook coefficient products.

### 7.3 Reassemble the Fp12 product

Compute H(A) and H(D), sharing \(B+C\) and \(B-C\), then set

$$
\begin{aligned}
A_{\rm new}&=2A+v^2H(D),\\
D_{\rm new}&=2D+vH(A).
\end{aligned}
$$

The result is exactly

$$
A_{\rm new}+D_{\rm new}w=2f\widehat L.
$$

There are two four-product calls:

$$
4M_2+4M_2=8M_2.
$$

Multiplication by v is cheap:

$$
v(t_0+t_1v+t_2v^2)=\xi t_2+t_0v+t_1v^2.
$$

### 7.4 The corresponding Python

This is the arithmetic in
[`_twice_mul_linear`](../../reference/normalized_miller_ref.py). The
excerpt assumes that `Fp6` is imported and all operands share the tower:

```python
def twice_mul_linear(t, B, C):
    m0 = t.c0 * B
    m3 = t.c2 * C
    s = t.c0 + t.c2
    mp = (s + t.c1) * (B + C)
    mm = (s - t.c1) * (B - C)
    return Fp6(
        (m0 + m3.mul_xi()).mul_int(2),
        mp - mm - m3.mul_int(2),
        mp + mm - m0.mul_int(2),
    )
```

The actual helper accepts precomputed `plus=B+C` and `minus=B-C` so the
two calls can share them. The surrounding `mul_line_scaled8` implements
the reassembly above.

For testing, the expected value is **twice** the dense product:

```python
got = nm.mul_line_scaled8(f, B, C)
check(got == scale12(f * dense_line, 2), "line identity")
```

Comparing `got` with the unscaled product would test the wrong contract.

## 8. Derive Fourier-based squaring

### 8.1 What the current general square does

For \(f=A+Dw\), \(w^2=v\),

$$
f^2=(A^2+vD^2)+2ADw.
$$

The current implementation obtains the real part from

$$
(A+D)(A+vD)-AD-vAD=A^2+vD^2.
$$

This uses two general Fp6 products:

$$
2M_6=12M_2=36M_p,
$$

plus additions and cheap v multiplications.

We now use the alternative view \(f=a+bw+cw^2\), with Fp4 coefficients.

### 8.2 First derive the polynomial square

Let

$$
F(z)=a+bz+cz^2,\qquad G(z)=F(z)^2.
$$

Write

$$
G(z)=q_0+q_1z+q_2z^2+q_3z^3+q_4z^4.
$$

The direct expansion is

$$
(q_0,q_1,q_2,q_3,q_4)
=(a^2,\ 2ab,\ b^2+2ac,\ 2bc,\ c^2).
$$

Our goal is to recover these coefficients, or a common multiple of them,
from five squarings.

### 8.3 Choose cheap evaluation points

Because \(u^2=-1\), the values \(1,-1,u,-u\) are fourth roots of unity.
They are distinct in this odd-characteristic field.

Compute

$$
\begin{aligned}
E_1&=(a+b+c)^2,\\
E_{-1}&=(a-b+c)^2,\\
E_u&=(a+ub-c)^2,\\
E_{-u}&=(a-ub-c)^2,\\
K&=4c^2.
\end{aligned}
$$

The fifth square supplies \(q_4=c^2\), the leading coefficient, or
“evaluation at infinity.” Four finite evaluations alone would not
determine an arbitrary degree-four polynomial.

This is a small, exact finite-field Fourier transform. It is not a
floating-point FFT.

### 8.4 Recover the coefficient combinations

Define

$$
\begin{aligned}
S&=E_1+E_{-1},\\
T&=E_u+E_{-u},\\
U&=E_1-E_{-1},\\
V&=-u(E_u-E_{-u}).
\end{aligned}
$$

Expanding each evaluation using \(u^2=-1\) gives

$$
\begin{aligned}
S&=2(q_0+q_2+q_4),\\
T&=2(q_0-q_2+q_4),\\
U&=2(q_1+q_3),\\
V&=2(q_1-q_3).
\end{aligned}
$$

Now add and subtract:

$$
\begin{aligned}
4q_0&=S+T-K,\\
4q_1&=U+V,\\
4q_2&=S-T,\\
4q_3&=U-V,\\
4q_4&=K.
\end{aligned}
$$

Again, we avoid division by keeping the common factor four.

### 8.5 Reduce to the three output coordinates

In the Fp4 tower, \(w^3=s\), so \(w^4=sw\).
Therefore

$$
G(w)=(q_0+sq_3)+(q_1+sq_4)w+q_2w^2.
$$

The scaled outputs are

$$
\boxed{
\begin{aligned}
h_0&=S+T-K+s(U-V),\\
h_1&=U+V+sK,\\
h_2&=S-T.
\end{aligned}}
$$

They satisfy the exact identity

$$
\boxed{h_0+h_1w+h_2w^2=4f^2.}
$$

No elliptic-curve assumption was needed for this local identity.
The pairing assumption is needed only when discarding the scale later.

### 8.6 Count the actual smaller squares

For \(a=a_0+a_1s\in\mathbb F_{p^4}\),

$$
a^2=(a_0^2+\xi a_1^2)
 +\left((a_0+a_1)^2-a_0^2-a_1^2\right)s.
$$

This implementation uses three S2, hence

$$
S_4=3S_2=6M_p.
$$

Our five evaluations cost

$$
5S_4=15S_2=30M_p.
$$

The multiplication count drops by \(6M_p\), or 16.7% of the original
square's multiplicative work. That is **not** a prediction of a 16.7%
runtime reduction: the linear work differs.

The Fourier implementation uses 19 external A4 operations, including
doublings, plus the cheap u and s operations. The five S4 calls have
their own internal additions.

### 8.7 Map the equation to the code

[`fft4_square`](scaled-square/gaussian_ref.py) uses these names:

| Mathematics | Python |
| --- | --- |
| \(E_1,E_{-1},E_u,E_{-u}\) | `pp, pm, pu, mu` |
| \(K=4c^2\) | `c4` |
| \(s\cdot a\) | `times_s(a)` |
| \(u\cdot a\) | `times_u(a)` |
| Fp4 addition/subtraction | `add(a,b)`, `sub(a,b)` |
| Fp4 squaring | `sqr4(a)` |

An Fp4 value is represented by a **pair of Fp2 values**. Ordinary Python
tuple `+` would concatenate tuples, so the explicit helpers matter.

The helper below implements the Fp4 square:

```python
def sqr4(a):
    x, y = a
    xx, yy = x.sqr(), y.sqr()
    xy2 = (x + y).sqr() - xx - yy
    return xx + yy.mul_xi(), xy2
```

After calculating `h0`, `h1`, and `h2`, reassemble the original storage
order:

```python
return Fp12(
    Fp6(h0[0], h2[0], h1[1]),
    Fp6(h1[0], h0[1], h2[1]),
)
```

Do not return `Fp6(h0[0], h1[0], h2[0])` as the first half; that would
confuse the Fp4 grouping with the original Fp6 grouping.

The C equivalent is `fft4_sqr` in
[kernels.h](scaled-square/kernels.h). `qsplit` and `qjoin` do the
regrouping. Temporaries make in-place calls such as `fft4_sqr(f,f)` safe;
the C tests check both in-place and separate-output use.

### 8.8 Why retain the other squaring candidates?

The archive also contains:

- `gaussian_square` / `gaussian_sqr`, using points
  \(\{0,1,-1,u,\infty\}\). It also costs 15S2 but uses 23 external A4
  operations and measured slower.
- `cubic_sqr`, a scaled Chung–Hasan-style control costing
  \(4S_4+M_4=33M_p\).

These distinguish a real arithmetic improvement from merely beating a
weak baseline. They also preserve the research trail.

## 9. Two hand-worked examples

These examples deliberately use smaller coefficients and a different field.
Work in

$$
\mathbb F_{13}[z]/(z^3-2).
$$

The nonzero cubes modulo 13 are \(1,5,8,12\), so 2 is not a cube and this
cubic polynomial is irreducible. For Fourier evaluation use \(u=5\),
since \(5^2\equiv-1\pmod{13}\).

The symbol u here is an ordinary element of F13. In the BLS12-381 tower,
u is the generator of Fp2. The interpolation identity needs only
\(u^2=-1\), so both settings are valid examples.

### 9.1 Four multiplications for a quadratic times a linear polynomial

Take

$$
T(z)=1+2z+3z^2,\qquad B+Cz=4+5z.
$$

The four products, reduced modulo 13, are

$$
m_0=4,\quad m_3=2,\quad m_+=2,\quad m_-=11.
$$

Substituting the reduction constant \(\xi=2\):

$$
\begin{aligned}
H_0&=2(4+2\cdot2)=3,\\
H_1&=2-11-2\cdot2=0,\\
H_2&=2+11-2\cdot4=5.
\end{aligned}
$$

Thus the scaled answer is \(3+5z^2\).

Check it by ordinary multiplication:

$$
T(z)(4+5z)=4+13z+22z^2+15z^3.
$$

After reducing \(z^3=2\) and all coefficients modulo 13, the product is
\(8+9z^2\). Twice that is \(3+5z^2\), as required.

### 9.2 Five squarings for the Fourier formula

Take the same \(F(z)=1+2z+3z^2\). The values are:

| Quantity | Value modulo 13 |
| --- | ---: |
| \(E_1=(1+2+3)^2\) | 10 |
| \(E_{-1}=(1-2+3)^2\) | 4 |
| \(E_u=(1+5\cdot2-3)^2\) | 12 |
| \(E_{-u}=(1-5\cdot2-3)^2\) | 1 |
| \(K=4\cdot3^2\) | 10 |

Then

$$
S=1,\quad T=0,\quad U=6,\quad V=10.
$$

The cubic reduction constant is now \(s=2\), so

$$
\begin{aligned}
h_0&=1+0-10+2(6-10)=9,\\
h_1&=6+10+2\cdot10=10,\\
h_2&=1-0=1.
\end{aligned}
$$

The scaled square is \(9+10z+z^2\).

Independently,

$$
F(z)^2=1+4z+10z^2+12z^3+9z^4.
$$

Reducing \(z^3=2,\ z^4=2z\) gives \(12+9z+10z^2\).
Multiplication by four gives \(9+10z+z^2\), matching the formula.

Both examples, their intermediate values, and an intentional corruption
control are executed by `toy_examples()` in
[tutorial_examples.py](tutorial_examples.py).

**Do not infer a pairing speedup from this toy field.** It illustrates
the algebra, not the cost model or final-exponent condition.

## 10. Combine the kernels safely

### 10.1 One line per iteration

Suppose the desired normalized update is

$$
f_{\rm next}=f^2\widehat L.
$$

Our composition returns

$$
\operatorname{scaled8}\bigl(\operatorname{fft4}(f),B,C\bigr)
=2(4f^2)\widehat L=8f^2\widehat L.
$$

It therefore has a scale of eight relative to the exact normalized update.

With k line applications after the square, it returns

$$
4\cdot2^k f^2\prod_{j=1}^k\widehat L_j.
$$

Relative to unnormalized original lines, there are also the factors
\((\ell_{a,j}y_{P_j})^{-1}\). They are nonzero Fp2 elements under the
validated-input preconditions.

### 10.2 Compare fairly with the existing 8M2 prototype

The existing normalized prototype already introduces the \(2^k\) factor.
When comparing the new loop with that prototype, the only additional
discrepancy is from the square:

$$
\alpha_0=1,\qquad
\alpha_{i+1}=4\alpha_i^2\pmod p.
$$

For example,

$$
\alpha_1=4,\quad \alpha_2=4^3,\quad \alpha_3=4^7,\quad
\alpha_n=4^{2^n-1}.
$$

The test helper `replay` tracks this scalar as a correctness diagnostic.
The timed C loop does **not** need to compute or store it.

After the final exponent, the two outputs are identical.

### 10.3 The e versus e-cubed convention

The definition-level routine computes \(e=f^E\).
The current fast BLS12 chain computes \(f^{3E}=e^3\), based on

$$
3\lambda=(x-1)^2(x+p)(x^2+p^2-1)+3,
\qquad \lambda=\frac{p^4-p^2+1}{r}.
$$

These are not identical raw target-group values in general. Cubing is
invertible on the order-r group here, so it preserves bilinearity,
nondegeneracy, and equality-to-one tests. Interoperability still requires
matching the same convention.

The discarded factors vanish under both E and 3E. The examples check
both conventions instead of confusing a convention mismatch with a kernel
failure. The BN fast convention in this repository is different: it
computes the exact exponent without the extra cube.

### 10.4 Leading operation counts

For a single prepared term in the current 64-square/69-line schedule:

| Path | Square cost | Line cost including evaluation | Online P setup | Total before final exponent |
| --- | ---: | ---: | ---: | --- |
| Original production arithmetic | \(36M_p\) | \(49M_p\) | none beyond the compared affine-input path | \(5685M_p\) |
| Normalized 8M2 prototype | \(36M_p\) | \(28M_p\) | \(1I_p+1M_p\) | \(4237M_p+1I_p\) |
| 8M2 plus Fourier square | \(30M_p\) | \(28M_p\) | \(1I_p+1M_p\) | \(3853M_p+1I_p\) |

The new square removes \(64\cdot6=384M_p\) per shared loop.
It does not remove any final-exponentiation operations.
These totals exclude additions, point validation, table preparation, and
other wrapper costs. Do not report their ratios as measured speedups.

## 11. Read the code in order

### 11.1 Suggested source map

| Read order | File and functions | Purpose |
| --- | --- | --- |
| 1 | [elips_ref.py](../../reference/elips_ref.py): `Fp2`, `Fp6`, `Fp12` | Textbook field values and operations |
| 2 | [pairing_ref.py](../../reference/pairing_ref.py): `line`, `miller`, `optimal_ate` | Independent dense pairing oracle |
| 3 | [normalized_miller_ref.py](../../reference/normalized_miller_ref.py): `prepare_g2`, `_twice_mul_linear`, `mul_line_scaled8` | First optimization |
| 4 | [gaussian_ref.py](scaled-square/gaussian_ref.py): `sqr4`, `fft4_square`, `replay` | Second optimization and shared replay |
| 5 | [tutorial_examples.py](tutorial_examples.py) | Hand examples, local identities, and complete pairing example |
| 6 | [kernels.h](scaled-square/kernels.h) | C implementation and comparison kernels |
| 7 | [experiment.c](scaled-square/experiment.c) | C correctness and timing harness |
| 8 | [inline_checks.py](inline_checks.py) | Exact symbolic identities and reconstructed exploratory checks |
| 9 | [analyze_bench.py](scaled-square/analyze_bench.py) | Summarize saved timing samples |

Production arithmetic is in [src/arith/fpx.c](../../../src/arith/fpx.c);
the production Miller loop is in
[src/pairing/miller.c](../../../src/pairing/miller.c). Neither was replaced
by these research kernels.

### 11.2 A complete pairing example

The following is the core flow in `pairing_example()` in the companion
script. Its imports and `check` helper are supplied by that script:

```python
P, Q, expected = nm._load_vectors(nm.DEFAULT_VECTORS)[0]
table = nm.prepare_g2(Q)

baseline = nm.miller_product_scaled8([(P, table)])
candidate, scale = replay([(P, table)], fft4_square)

check(candidate == scale12(baseline, scale), "whole-loop scale")
e = nm.final_exponentiate(candidate)
check(e == expected, "known-answer vector")
```

The first comparison checks the exact raw discrepancy. The second checks
the final value against an independently generated committed vector.

The loader starts with `_` because it is a private test utility, not a
stable application API. `replay` also assumes trusted valid inputs; it
does not implement the validation policy of a production pairing API.

### 11.3 Reading C safely

- `fp_t`, `fp2_t`, and related types are arrays. Pass them directly to the
  arithmetic routines, not as pointer-to-pointer values.
- Use conversion helpers and `fp_set_one`; respect Montgomery form.
- Check both aliasing cases: `result != input` and `result == input`.
- Keep the existing exact-raw API separate from scaled research paths.
- Do not copy only the BLS12 replay into a BN implementation. BN needs
  its additional Frobenius correction lines.

## 12. Run and reproduce the experiments

All commands below assume the **repository root** as the current directory.
No step requires committing, pushing, installing the library system-wide,
or deleting an old experiment.

### 12.1 First run: no C compiler and no extra Python packages

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  tools/research/pairing/tutorial_examples.py --skip-pairing

PYTHONDONTWRITEBYTECODE=1 python3 \
  tools/research/pairing/tutorial_examples.py
```

The first command checks the short examples and full-size field identities.
The second also checks a complete pairing KAT, the raw scale relation,
the e-cubed convention, and cancellation.

Successful output ends with:

```text
All requested tutorial examples passed.
```

These checks explicitly raise on failure and still work with `python3 -O`.

### 12.2 Broader Python checks

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/reference/normalized_miller_ref.py
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/scaled-square/gaussian_ref.py
```

These include field edge cases, known-answer pairing vectors, shared
products, operation counters, and negative controls.

Optional exact symbolic checks require SymPy:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/research/pairing/inline_checks.py
```

If SymPy is missing, install it into a separate environment of your choice,
not by changing system Python to make this tutorial work. The standard-library
examples above do not require it.

### 12.3 Build the research programs

Requirements: C11 compiler, CMake, and GMP development files.
Choose a new build directory if the example path is already in use:

```sh
cmake -S tools/research/pairing -B /tmp/elips-note-release \
  -DCMAKE_BUILD_TYPE=Release -DELIPS_BUILD_TESTS=OFF -DELIPS_WERROR=ON
cmake --build /tmp/elips-note-release -j 4
ctest --test-dir /tmp/elips-note-release --output-on-failure
```

`ELIPS_BUILD_TESTS=OFF` disables the main library's test targets in this
standalone build. It does **not** disable the two research CTest tests.
This is not a command to run the complete repository suite.

The executables are:

```text
normalized_line_experiment
scaled_square_experiment
primitive_ratios
```

For memory and undefined-behavior checks:

```sh
cmake -S tools/research/pairing -B /tmp/elips-note-asan \
  -DCMAKE_BUILD_TYPE=Asan -DELIPS_BUILD_TESTS=OFF -DELIPS_WERROR=ON
cmake --build /tmp/elips-note-asan -j 4
ctest --test-dir /tmp/elips-note-asan --output-on-failure
```

Do not compare sanitizer timings with optimized Release timings.

### 12.4 Benchmark without overwriting evidence

```sh
/tmp/elips-note-release/scaled_square_experiment --bench
/tmp/elips-note-release/normalized_line_experiment --bench
PYTHONDONTWRITEBYTECODE=1 python3 \
  tools/research/pairing/scaled-square/analyze_bench.py
```

The last command reads the preserved timing logs. It does not rerun the
benchmark. If saving a new run, give it a new file name and retain previous
results. Avoid running builds, tests, or another benchmark concurrently.

The archive index, [README.md](README.md), also lists the optional legacy
norm and lattice experiments. The lattice script defaults to quick unit
tests; `--search` enables a potentially slow heuristic computation.

## 13. Interpret the measurements

### 13.1 Keep the two comparisons separate

The earlier normalized-line experiment compared the original 15M2 sparse
kernel with the 8M2 candidate. It reported roughly 20–31% less prepared
Miller time across its tested 1-, 2-, and 8-term workloads.
Its original raw logs and conditions are in
[normalized-line/RESEARCH.md](normalized-line/RESEARCH.md).

The later squaring experiment held the **8M2 line kernel fixed**.
Its gains are additional to that baseline, not another comparison with
the original 15M2 path.

| Workload | Run 1: less median time | Run 2: less median time |
| --- | ---: | ---: |
| General squaring | 10.15% | 10.18% |
| Prepared Miller, 1 term | 5.71% | 5.37% |
| Prepared Miller + final exponent, 1 term | 1.88% | 1.61% |
| Prepared Miller, 2 terms | 3.24% | 3.66% |
| Prepared Miller + final exponent, 2 terms | 2.05% | 1.19% |
| Prepared Miller, 8 terms | 0.97% | 0.79% |
| Prepared Miller + final exponent, 8 terms | 1.07% | 0.62% |

The platform was Apple M2 Max, Apple Clang 21, GMP 6.3.0, with an
optimized Release build. Each run used 31 samples per method/workload,
rotating method order and warming up before measurement.

See [RESULTS.md](scaled-square/RESULTS.md) for the full method, comparisons
against the asymmetric cubic control, and exclusions. In particular, the
small whole-pairing improvements over that stronger control were noisy.

### 13.2 Why the whole pairing improves less

If a fraction q of a workload becomes a times as fast, the idealized
new/old time ratio is

$$
(1-q)+\frac{q}{a}.
$$

This is Amdahl's law. Accelerating the Miller squarings does not accelerate
the unchanged final exponentiation, point validation, or preparation.

Multi-pairing shares one square across many lines. Saving a fixed amount
per square therefore becomes a smaller percentage as the number of terms
increases.

### 13.3 What a careful benchmark must say

Always identify:

1. The exact baseline and candidate.
2. Which setup and validation costs are included.
3. Whether final exponentiation is included.
4. The number of terms and reuse assumptions.
5. Compiler, build flags, processor, and measurement method.
6. Variation and adverse results, not only the best sample.

The reported measurements exclude Q preparation, subgroup checks,
point-to-affine conversion, and the public API's chunk wrapper. Online P
normalization is included. These are **prepared arithmetic-path**
measurements, not end-to-end BLS signature verification measurements.

## 14. Correctness and security checklist

### 14.1 Layer the evidence

| Layer | Question | Evidence in the archive |
| --- | --- | --- |
| Polynomial identity | Are the equations correct for arbitrary coefficients? | Exact symbolic expansion and reduction |
| Field implementation | Does code implement that identity? | Random/edge field comparisons against dense operations |
| Loop composition | Are scale factors handled across every iteration? | Exact raw scale recurrence |
| Pairing value | Is the final answer still the intended pairing? | Committed independent KATs and exact final exponent |
| Multi-pairing | Does sharing the accumulator preserve products? | Empty, multi-term, and cancellation tests |
| Memory behavior | Are array accesses and aliasing safe in tested paths? | In-place checks, ASan, UBSan |
| Performance | Is the complete compared workload faster? | Repeated interleaved C measurements |

The Fourier reference checks cover 1,009 field cases on each of three
curves, all four BLS12-381 pairing vectors, and shared products.
The C comparison covers 4,000 field cases and multiple complete
BLS12-381 pairing workloads for all four square methods.

Those counts describe the existing experiment suites. The new companion
tutorial is smaller and designed to explain the tests, not replace them.

### 14.2 Common mistakes

- Comparing `4*f^2` with `f^2` and declaring the formula incorrect.
- Removing a factor from only one coordinate instead of all coordinates.
- Dropping a zero denominator or silently using `inv(0)`.
- Using subgroup-specific squaring on a general Miller accumulator.
- Assuming Python tuples have extension-field arithmetic.
- Confusing the C storage order with increasing powers of w.
- Treating two values that agree only after final exponentiation as
  interchangeable in an exact-raw API.
- Treating a low operation count as a measured speedup.
- Treating a failed heuristic search as an impossibility proof.

### 14.3 Security limits

The Python prototypes are variable-time and use deterministic test inputs.
The research C harness is also not a reviewed application boundary.
No constant-time leakage measurement or production security review was
performed for these new paths.

A production integration would need, at minimum:

- input and subgroup validation, including malformed encodings;
- explicit identity and zero-denominator behavior;
- a separate or versioned contract for changed raw values;
- constant-time review of the chosen implementation;
- broader tests, cross-implementation checks, and target-platform benchmarks.

ASan and UBSan do not establish cryptographic security or absence of
timing leakage. Mathematical identities do not establish memory safety.

## 15. Other research paths

These are preserved to show how to document promising and unsuccessful
directions without overstating them.

### 15.1 Tangent coefficients from the curve equation

With homogeneous twist coordinates \(x=X/Z,\ y=Y/Z\), the curve equation is

$$
Y^2Z=X^3+b_{\rm twist}Z^3.
$$

The tangent coefficients used by the compared implementation are

$$
\ell_a=2\xi YZ^2,\quad
\ell_b=-3X^2Z,\quad
\ell_c=3X^3-2Y^2Z.
$$

Using the curve equation,

$$
(\ell_a,\ell_b,\ell_c)
=Z\left(2\xi YZ,\ -3X^2,\ Y^2-3b_{\rm twist}Z^2\right).
$$

For a finite nonidentity state, the common nonzero Z factor can be
discarded in a reduced-pairing path. Some remaining expressions are
already computed by point doubling, offering reuse.

The original experiment checked this identity on 32 projectively rescaled
points. It did not implement and benchmark a complete optimized doubling
routine. This reduction is separate from the prepared-loop speedups above.

### 15.2 Prescribed norms

The norm from \(\mathbb F_{p^{12}}\) to \(\mathbb F_{p^6}\) is

$$
N_{12/6}(f)=f\,f^{p^6}.
$$

To \(\mathbb F_{p^4}\), it is

$$
N_{12/4}(f)=f\,f^{p^4}f^{p^8}.
$$

The denominator-retaining experiment tests, for \(T=[m]Q\), identities of
the form

$$
N_{12/6}(f)=\frac{(x_P-x_Q)^m}{x_P-x_T},\qquad
N_{12/4}(f)=\frac{(y_P-y_Q)^m}{y_P-y_T}.
$$

These use untwisted coordinates and a particular normalized, complete
Miller function. Do not substitute an arbitrary denominator-eliminated
or scaled accumulator and assume the same right-hand sides still apply.

The idea was to exploit these constraints to represent less state.
The script verified small m cases on BLS12-381 and BN-462; it did not
produce a cheaper compressed recurrence or a cheap runtime self-check.
Compressed pairings and torus representations also have substantial
published prior art.

### 15.3 Final-exponent searches

The hard exponent is

$$
\lambda=\frac{p^4-p^2+1}{r}.
$$

The preserved lattice experiment seeks a small-coefficient representation

$$
3\lambda\equiv\sum_{i,j} c_{ij}|x|^i p^j
\pmod{p^4-p^2+1}.
$$

The modulus is the cyclotomic subgroup order because that is where the
easy-part output lies.

LLL reduces a lattice basis; Babai's method returns an approximate
closest-vector candidate. Recovering a known small representation is a
useful positive control. Returning large coefficients for another basis
does not prove that every representation has large coefficients.

Even a short representation is not automatically a fast exponentiation
program. Common subexpressions, operation order, inversions, and
decompression costs all matter.

## 16. Exercises and next experiments

### Exercise 1: recover the missing factor

The line kernel returns \(2f\widehat L\); the square returns \(4f^2\).
What is the output with two line applications after one square?

**Answer:** \(16f^2\widehat L_1\widehat L_2\).

### Exercise 2: spot an invalid comparison

If the optimized raw Miller output differs from the baseline, must the
pairing be wrong?

**Answer:** no. Check the promised scale relation and compare the final
exponentiated values. Conversely, raw inequality alone says nothing about
whether the discrepancy is permitted.

### Exercise 3: explain the fifth evaluation

Why do four fourth roots of unity not suffice to interpolate \(F(z)^2\)?

**Answer:** it has degree up to four, hence five coefficients. At fourth
roots of unity, \(z^4=1\), so the constant and degree-four terms alias.
The independent value \(c^2=q_4\) separates them.

### Exercise 4: test your own change

Alter the sign of V in a **new copy** of the Fourier experiment.
Check that the toy example and full-size field comparisons reject it.
Keep the failed experiment and label it, rather than overwriting a
working version.

### Exercise 5: investigate a performance hypothesis

Can fewer copies or better reuse of additions reduce the measured cost
without changing the polynomial identity?

Before implementing, write down:

1. the exact output contract;
2. the old and new operation counts, including additions;
3. an independent correctness check;
4. the benchmark that would disprove the speed hypothesis.

Useful future directions include combined square-and-line circuits,
carefully bounded lazy reduction, preparation amortization, and a
different Miller digit schedule. These are research opportunities,
not promises of improvement or novelty.

## 17. Further reading

Read these for background and attribution, not as evidence that our exact
implementation is novel:

1. **Miller's algorithm, introductory notes.**
   [Stanford PBC notes](https://crypto.stanford.edu/pbc/notes/ep/miller.html).
   Explains line functions and the basic loop.
2. **Devegili, Ó hÉigeartaigh, Scott, and Dahab,
   Multiplication and Squaring on Pairing-Friendly Fields.**
   [ePrint 2006/471](https://eprint.iacr.org/2006/471).
   Directly relevant: its section 2 explicitly discards subfield scale
   factors to clear Toom-Cook denominators. This principle predates this work.
3. **Chung and Hasan, Asymmetric Squaring Formulae.**
   [Paper](https://www.lirmm.fr/arith18/papers/Chung-Squaring.pdf).
   Background for the stronger cubic-squaring control.
4. **Costello and Stebila, Fixed Argument Pairings.**
   [ePrint 2010/342](https://eprint.iacr.org/2010/342).
   Background for prepared pairing workloads.
5. **Naehrig, Barreto, and Schwabe,
   On compressible pairings and their computation.**
   [Paper](https://cryptojedi.org/papers/ocpatc.pdf).
   Important prior art for norm/torus representations during pairing computation.
6. **Granger and Scott,
   Faster squaring in the cyclotomic subgroup of sixth degree extensions.**
   PKC 2010. The production final exponentiation uses this subgroup-specific
   idea; it is distinct from the general scaled square derived here.

For the exact local results, read
[scaled-square/RESULTS.md](scaled-square/RESULTS.md), the saved logs, and
[CORRECTIONS.md](CORRECTIONS.md) together.

The defensible research conclusion is:

> A final computation's invariances can expose cheaper intermediate
> arithmetic. Here that leads to explicit, tested formulas and modest
> incremental speedups on one implementation. It does not prove an optimal
> pairing algorithm, a universally faster implementation, or a novel
> cryptographic result.
