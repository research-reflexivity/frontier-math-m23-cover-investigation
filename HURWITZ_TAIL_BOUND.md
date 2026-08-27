# The uniform Hurwitz tail bound

## Status

The analytic comparison has three logically separate parts.

1. `verification/verify_hurwitz_tail_geometry.py` proves, with Arb balls, that
   every point of the two half-triangles is within pseudohyperbolic distance
   `0.471` of one of the 47 atlas centres.  The computed upper endpoint is
   `0.4708683714328407`.
2. `scripts/analyze_hurwitz_tail_stability.py` certifies a normalized left
   inverse for the modes `0,...,60`, bounds the coupling to the remaining
   modes by Parseval and a Schur complement, and closes the outer-sup-norm
   bootstrap at radius `R=0.99`.
3. The final model computation must bound **all** `Q=1280` DFT output modes.
   Checking only modes through the polynomial cutoff is insufficient.  The
   degree-480 models have residuals below `8.5e-93` in modes `0,...,480`, but
   their full 1280-mode residual norms are about `2e-67`.  Those high residuals
   are the honest degree-480 truncation, not interval noise.  A longer model is
   therefore used for the final digit count.

The first two parts pass for all seven Hurwitz classes.  This note records the
inequalities so that the third part can be substituted without changing the
proof.

## Scaled Taylor coefficients

On the chart with centre `z_j`, write

```text
g_j(w) = f(z_j(w))/(1-w)^2 = sum_{n>=0} b_{j,n} w^n,
c_{j,n} = b_{j,n} rho^n,
rho = 0.72.
```

Thus the polynomial evaluated by the code is

```text
sum_n c_{j,n} (w/rho)^n.
```

Put

```text
M_R(f) = max_j sup_{|w|<=R} |g_j(w)|,       R = 0.99.
```

Cauchy's inequality gives, with `r=rho/R`,

```text
|c_{j,n}| <= M_R(f) r^n.                    (1)
```

This is the only Taylor-coefficient estimate used below.

## Why the 47 charts cover every point

The fundamental triangle is split into two half-triangles and placed in the
Klein disk.  A Euclidean triangle in the Klein disk is hyperbolically
geodesically convex.  For every rational mesh cell the verifier chooses a
mesh vertex `v`, proves that all three vertices are in a hyperbolic ball of
pseudohyperbolic radius `epsilon` about `v`, and hence proves the same for the
whole cell.  If `v` is at pseudohyperbolic distance `d` from the selected
atlas centre, the pseudohyperbolic triangle inequality gives

```text
delta(cell, centre) <= (d+epsilon)/(1+d epsilon).
```

The largest Arb upper endpoint over the 51,200 cells is
`0.4708683714328407`, so the round value

```text
delta = 0.471 < rho
```

is valid uniformly, not merely on the mesh vertices.

For example, the pointwise target Taylor tail after mode `N` is therefore

```text
M_R(f) (delta/R)^(N+1)/(1-delta/R).          (2)
```

At `N=480` the coefficient multiplying `M_R(f)` in (2) is less than
`1.268e-155`.

## The finite collocation operator

Let `q=0,...,Q-1` index the source points `w_s=rho*zeta_Q^q`.  Triangle
reduction chooses a target chart `t(i,q)`, a target coordinate `w_t`, and a
weight-two transition factor `F(i,q)`.  Set `u(i,q)=w_t/rho`.  The exact
sampled equation is

```text
sum_n c_{i,n} zeta_Q^(qn)
  = F(i,q) sum_n c_{t(i,q),n} u(i,q)^n.      (3)
```

Applying the `1/Q`-normalized DFT turns (3) into a coefficient equation.  For
indices below `Q`, the source side is the identity.  Coefficients with indices
at least `Q` appear as source aliases.

Take `L=60` and split the coefficients below `Q` into

```text
x: modes 0,...,L,
y: modes L+1,...,Q-1.
```

The four anchor equations are appended to the low equations.  If `A` denotes
that augmented low-low block, the stored binary64 inverse is used only as a
candidate.  Its residual, matrix-product roundoff, exact roots of unity, and
the difference between binary64 routes and Acb routes are all included in an
a posteriori error `eta<1`.  Consequently

```text
sigma_min(A) >= (1-eta)/||B||.               (4)
```

Across the seven classes the lower bounds in (4) lie between `6.60e-4` and
`9.01e-4`.

## Parseval bounds and the Schur complement

For a set `E` of input modes, Parseval gives the routewise bound

```text
||T_E||^2 <= max_j (1/Q) sum_{t(i,q)=j}
  |F(i,q)|^2 sum_{n in E} |u(i,q)|^(2n).     (5)
```

The single factor `1/Q` in (5) is important: the DFT is normalized by `1/Q`,
and Parseval is being applied to the squared norm.  Acb-enclosed route radii
and factors give

```text
l = ||T_{0,...,L}|| <= 1.38e3,
h = ||T_{L+1,...}|| <= 2.35e-11.
```

Writing the high-high block as `T_hh-I`, its inverse norm is at most
`1/(1-h)`.  Eliminating `y` changes the low singular gap by at most

```text
h*l/(1-h).
```

All seven certified Schur margins

```text
m = sigma_min(A) - h*l/(1-h)                 (6)
```

are positive, again between `6.60e-4` and `9.01e-4`.  If a full residual is
split as `(r_low,r_high)`, (6) gives the sharper low-coefficient estimate

```text
||x|| <= (||r_low|| + h ||r_high||/(1-h))/m. (7)
```

In particular, a residual supported only in high output modes reaches the
leading coefficients with factor at most

```text
h/((1-h)m) < 2.79e-8.                        (8)
```

This is much sharper than multiplying the full residual by the norm of the
whole inverse.

## Closing the outer-sup-norm bootstrap

For a form normalized by one of the four anchor vectors, (1) gives

```text
C_tail = sqrt(47) r^(L+1)/sqrt(1-r^2) < 3.657e-8
```

for the global `l2` norm of the coefficients beyond `L`.  The low source
aliases have norm at most

```text
C_alias = sqrt(47*(L+1)) r^Q/(1-r^Q) < 5.03e-176.
```

Evaluation at an inner point of radius at most `delta` has low-mode norm

```text
E_low = sqrt(sum_{n=0}^L (delta/rho)^(2n)) < 1.323.
```

The weight-two transition formula is

```text
|F| = (y_t/y_s) (1-|w_t|^2)/(1-|w_s|^2).
```

Arb gives a uniform outer-chart bound below `79628`; the proof uses `80000`.
Combining this with (4), (5), and the pointwise tail (2) at `N=L` gives

```text
M_R <= 80000 * [ E_low/sigma_min(A)
                 * (1 + (h*C_tail+C_alias) M_R)
                 + (delta/R)^(L+1)/(1-delta/R) M_R ].
```

The coefficient of `M_R` on the right is below `1.1e-10` in every class, so
the inequality is contractive.  The resulting anchor-normalized upper bounds
are below `1.7e8`.

There is no circular assumption that the four anchors normalize the exact
space.  If an exact differential has all four anchors equal to zero, the
constant `1` disappears from the displayed inequality.  Contractivity then
forces `M_R=0`, so the differential is zero.  The anchor map is therefore
injective on the known four-dimensional space of holomorphic differentials,
and hence is an isomorphism.  The four exact anchor-normalized forms used in
the comparison consequently exist and are unique.

The leading four-by-four branch-jet matrix is also certified invertible.  A
small Acb ball matrix is checked by its own a posteriori inverse, and the
low-coefficient error is subtracted from its singular gap.  The branch-
normalized bounds are below `3.0e10`.

## External forcing at Q=1280

The square `Q`-mode system treats modes `0,...,Q-1` as unknowns.  Hence its
genuine target tail starts at `Q`, not at the polynomial degree used to
propose the solution.  With `s=max|w_t|/R`, Parseval and (1) give

```text
target external <= sqrt(47) max|F| s^Q/(1-s) M_R,
source aliases  <= sqrt(47) r^Q/
                   (sqrt(1-r^2) (1-r^Q)) M_R.       (9)
```

Here `max|w_t|<0.459871`, so the first term in (9) is far below binary64's
normal range.  The source term dominates; after anchor normalization the
combined external forcing is below `2e-168` (the implementation retains a
larger positive number instead of underflowing the target term to zero).

## What the degree-480 computation does and does not prove

The degree-480 Acb refinement makes the residual in output modes
`0,...,480` about `1e-93`.  It does not solve output modes `481,...,1279`.
Evaluating all of them gives the following full residual scale:

```text
class  maximum full 1280-mode residual
1      1.95e-67
2      2.06e-67
3      1.68e-67
4      2.05e-67
5      2.01e-67
6      2.42e-67
7      1.71e-67
```

Equation (8) suppresses their effect on the leading coefficients by about
eight further decimal orders, but the subsequent branch normalization loses
some of that gain.  Thus the degree-480 batch remains an excellent independent
numerical cross-check, but it is not the cleanest batch on which to base a
high-digit uniform certificate.

The final computation therefore extends the solved polynomial cutoff to
`N=700` while keeping `Q=1280`.  Its outward-rounded results are:

```text
class  maximum all-mode residual  maximum first-20 jet-row error  digits
1      3.982e-93                  1.450e-76                       75
2      3.244e-93                  6.535e-77                       76
3      2.661e-93                  7.047e-77                       76
4      6.755e-93                  8.766e-77                       76
5      3.806e-93                  4.710e-77                       76
6      3.471e-93                  3.539e-77                       76
7      5.204e-93                  2.345e-76                       75
```

For completeness, the last propagation step is as follows.  Let `E` be the
Frobenius bound obtained by applying (7) to the four columns, let `L` be the
exact anchor-normalized leading branch matrix, and let `Lhat` be the computed
one.  Then

```text
||Lhat-L|| <= Delta = E/rho^3.
```

The separately certified lower bound for `sigma_min(L)` proves

```text
||Lhat^(-1)|| <= 1/(sigma_min(L)-Delta).
```

If `G` and `Ghat` are the exact and computed branch-normalized series, the
identity

```text
Ghat-G = (D-G*(Lhat-L))*Lhat^(-1)
```

combines the raw coefficient error `D` with the Cauchy estimate
`||G_n||_2 <= 2 M_R/R^n`.  Applying it for `n=0,...,19` gives the final column
of the table.  Thus every Taylor row used by the canonical quadric and Petri
cubic reconstruction has at least 75 certified decimal digits.
