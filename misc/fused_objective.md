# The fused normal-equations objective

A technical account of what `beta_method = "crossprod"` computes, why the fused
implementation is algebraically identical to the two-pass one it replaced, and
exactly where and by how much the two differ in floating point.

Every numerical claim below is measured on the SPX calibration
($n = 2428$, $t_c = 2717$, $m = 0.9$, $\omega = 6.1597$, R 4.5.2, Apple Silicon,
unit roundoff $u = 2^{-53} \approx 1.11\times10^{-16}$). The script is
`B1_fused_math.R`.

## 0. What is actually being computed

One number: the sum of squared errors.

Everything below — the design matrix, the basis, the projection, the normal
equations — is scaffolding for evaluating a single scalar

$$S(t_c, m, \omega) \;=\; \sum_{i=1}^{n}\big(y_i - \hat y_i\big)^2$$

as fast and as accurately as it can be evaluated. That scalar is the objective
`optim()` minimises over the three nonlinear parameters, and it is the *only*
thing the fused routine returns; the coefficients $\hat\beta$ are computed inside
it and discarded.

The reason any of this is worth the trouble is the call count. A single SPX F2
calibration performs **27,402** linear solves, of which **26,890** are objective
evaluations — 500 candidate $t_c$ values, each driving an L-BFGS-B search whose
finite-difference gradients cost roughly 54 objective evaluations apiece. The
remaining 512 recover the coefficients once per candidate after its optimisation
has converged, and are negligible. So a saving of a few microseconds per
evaluation is a saving of seconds per fit, and the numerical properties of this
one expression are the numerical properties of the whole calibration.

---

## 1. Separable structure

The Filimonov parameterisation (mode 0) of the LPPLS model is

$$f(t) = A + \tau^{m}\big[B + C_1\cos(\omega\log\tau) + C_2\sin(\omega\log\tau)\big],
\qquad \tau = t_c - t > 0 .$$

Split the parameters into $\beta = (A, B, C_1, C_2)' \in \mathbb{R}^4$ and
$\theta = (t_c, m, \omega)$. Distributing $\tau^m$ gives

$$f(t) = A\cdot 1 \;+\; B\cdot\tau^{m} \;+\; C_1\cdot\tau^{m}\cos(\omega\log\tau)
\;+\; C_2\cdot\tau^{m}\sin(\omega\log\tau),$$

so $f$ is **linear in $\beta$ for fixed $\theta$**. Collecting the four
coefficient functions at the sample times $t_1 < \cdots < t_n < t_c$ gives the
design matrix

$$X(\theta) = \big[\,\mathbf{1},\; \boldsymbol\tau^{m},\;
\boldsymbol\tau^{m}\!\odot\cos(\omega\log\boldsymbol\tau),\;
\boldsymbol\tau^{m}\!\odot\sin(\omega\log\boldsymbol\tau)\,\big]
\in \mathbb{R}^{n\times 4},$$

with $\odot$ elementwise, so that the vector of model values is exactly
$X(\theta)\beta$.

This is what the calibration exploits: rather than optimising over all seven
parameters, it minimises the **concentrated** (profiled) objective

$$S(\theta) \;=\; \min_{\beta\in\mathbb{R}^4}\ \lVert y - X(\theta)\beta\rVert_2^2$$

over $\theta$ alone — the variable-projection idea of Golub and Pereyra. The
inner minimisation is ordinary least squares and is solved in closed form, which
is why the solver runs once per objective evaluation and dominates the cost of a
fit.

---

## 2. The basis

Write $\phi_1 = 1$, $\phi_2 = \tau^m$, $\phi_3 = \tau^m\cos(\omega\log\tau)$,
$\phi_4 = \tau^m\sin(\omega\log\tau)$ as functions on $(0,\infty)$, and
$\mathcal V = \operatorname{span}\{\phi_1,\dots,\phi_4\}$.

**Proposition 1.** *If $m \neq 0$ and $\omega \neq 0$, the functions
$\phi_1,\dots,\phi_4$ are linearly independent on any subinterval of
$(0,\infty)$ with nonempty interior.*

*Proof.* Suppose $\sum_j a_j\phi_j$ vanishes on such a window $W$. Each $\phi_j$
is real-analytic on $(0,\infty)$ (as $\tau^m = e^{m\log\tau}$, composed with
$\cos$ and $\sin$), so by the identity theorem the combination vanishes on all of
$(0,\infty)$. Substituting $u = \log\tau \in \mathbb{R}$ and writing
$a_3\cos(\omega u) + a_4\sin(\omega u) = R\cos(\omega u - \varphi)$ with
$R = \sqrt{a_3^2 + a_4^2}$,

$$a_1 + e^{mu}\big(a_2 + R\cos(\omega u - \varphi)\big) = 0 \qquad \forall u\in\mathbb{R}.$$

Take $m > 0$ and let $u \to -\infty$: $e^{mu}\to 0$ while the bracket is bounded
by $|a_2| + R$, so the second term vanishes and $a_1 = 0$. (For $m < 0$, take
$u \to +\infty$.) With $a_1 = 0$ and $e^{mu} > 0$, we get
$a_2 + R\cos(\omega u - \varphi) \equiv 0$; since $\omega\neq0$ the cosine attains
both $+1$ and $-1$, forcing $R = 0$ and then $a_2 = 0$. Finally $R = 0$ gives
$a_3 = a_4 = 0$. $\blacksquare$

Both hypotheses are sharp, and the failure modes are the obvious ones:

| | collision | $\operatorname{rank} X$ |
|---|---|---|
| $m = 0$ | $\phi_2 = \phi_1$ | 3 |
| $\omega = 0$ | $\phi_3 = \phi_2$, $\phi_4 \equiv 0$ | 2 |

Both lie outside the admissible box ($0.1 \le m \le 0.9$, $6 \le \omega \le 13$),
where $\operatorname{rank} X = 4$ throughout.

**Proposition 2.** *Under Proposition 1, $\operatorname{rank} X = 4$ provided the
sample contains four points $s_1,\dots,s_4$ with $\det[\phi_j(s_i)] \neq 0$.*

*Proof.* Let $E:\mathcal V\to\mathbb{R}^n$ be evaluation at the sample points, so
$\operatorname{col}(X) = E(\mathcal V)$ and
$\operatorname{rank}X = 4 - \dim\ker E$. A nonzero $f\in\mathcal V$ is
real-analytic and not identically zero, so its zeros are isolated and finite in
the window; hence $\ker E = \{0\}$ unless the sample lies inside that finite zero
set. Proposition 1 makes the evaluation map on the whole window rank 4, so such a
quadruple exists. $\blacksquare$

This is not marginal in the combinatorial sense. Over the SPX window
$u=\log\tau$ spans length $2.240$, so the oscillatory factor completes only
$L\omega/2\pi = 2.20$ periods; 2,000 random nonzero elements of $\mathcal V$ had
at most **6** zeros (median 4), and four evenly spaced sample points already give
rank 4. With $n = 2428$ the columns are a basis with enormous margin
*combinatorially*. The margin that is not large is numerical, and $\kappa$ is what
measures it (§6).

---

## 3. The objective is a projection

With $X$ of full column rank, $S(\beta) = (y - X\beta)'(y - X\beta)$ expands to
$y'y - 2\beta'X'y + \beta'X'X\beta$, so

$$\nabla_\beta S = -2X'y + 2X'X\beta, \qquad \nabla^2_\beta S = 2X'X \succ 0 .$$

The stationary point is therefore the unique minimiser and satisfies the **normal
equations**

$$X'X\,\hat\beta = X'y, \qquad\text{so}\qquad \hat\beta = (X'X)^{-1}X'y,$$

and the concentrated objective is the squared distance from $y$ to the
four-dimensional subspace $\operatorname{col}(X)$:

$$S(\theta) = \lVert y - X\hat\beta\rVert^2 = \lVert (I - P)y\rVert^2, \qquad
P = X(X'X)^{-1}X' .$$

Verified numerically: $\lVert y - X\hat\beta\rVert^2$ and $\lVert (I-P)y\rVert^2$
agree to all printed digits ($5.2867886392$), and the residual is orthogonal to
every column, $\max_j |(X'r)_j| = 1.8\times10^{-9}$ against column norms of order
$10^{4}$.

---

## 4. The two implementations

Both compute the same $S$. They differ only in how the fitted values
$\hat y = X\hat\beta$ are formed.

**Algorithm A (two-pass, the previous implementation).**

1. Build $X$: compute $\boldsymbol\tau$, $\boldsymbol\tau^m$, $\log\boldsymbol\tau$, $\cos$, $\sin$.
2. Solve $X'X\hat\beta = X'y$.
3. *Discard $X$.* Recompute $\boldsymbol\tau$, $\boldsymbol\tau^m$, $\log\boldsymbol\tau$, $\cos$, $\sin$ and evaluate the closed form
   $$\hat y_i^{A} = \hat A + \tau_i^{m}\big(\hat B + \hat C_1 c_i + \hat C_2 s_i\big).$$
4. $S = \sum_i (y_i - \hat y_i^{A})^2$.

**Algorithm B (fused, current).**

1. Build $X$ once.
2. Solve $X'X\hat\beta = X'y$.
3. $\hat y^{B} = X\hat\beta$, i.e.
   $$\hat y_i^{B} = \hat A\cdot 1 + \hat B f_i + \hat C_1 g_i + \hat C_2 h_i,
     \qquad f_i = \tau_i^m,\ g_i = f_i c_i,\ h_i = f_i s_i .$$
4. $S = \sum_i (y_i - \hat y_i^{B})^2$.

**Proposition 3.** *In exact arithmetic $\hat y^{A} = \hat y^{B}$, hence the two
algorithms return the same $S$.*

*Proof.* Immediate from distributivity:
$\tau^m(B + C_1c + C_2s) = B\tau^m + C_1\tau^mc + C_2\tau^ms$. $\blacksquare$

The saving is that steps 1 and 3 of Algorithm A evaluate the same transcendental
quantities twice. Algorithm B walks the series once.

Note what step 4 returns in both cases: the scalar $S$, and nothing else. The
coefficients $\hat\beta$ are a means to it and are dropped on return. They are
recovered separately, once per candidate $t_c$ after its optimisation converges,
by `create_beta_calculator()` — which is why the package keeps the two factories
apart, and why the redundant re-solve there (512 of the 27,402, about 39 ms per
SPX fit) is not worth eliminating.

### What is shared bitwise

This is the important structural point. Both algorithms build $X$ with the same
expression in the same order, so $X$, $X'X$, $X'y$ and the LAPACK solve are
identical operations on identical inputs. Measured:

$$\max_j \big|\hat\beta^A_j - \hat\beta^B_j\big| = 0 \quad\text{(bitwise identical)} .$$

**The fusion does not change the solve at all.** It changes only the evaluation
of $\hat y$ from $\hat\beta$. Any statement about conditioning, solver accuracy or
the choice among `beta_method` values is therefore untouched by it (§6).

---

## 5. Floating-point analysis of the difference

Under the standard model $\mathrm{fl}(a \circ b) = (a\circ b)(1+\delta)$,
$|\delta| \le u$, the two groupings are

$$\hat y_i^{A}:\quad \hat A \oplus \big(f_i \otimes (\hat B \oplus (\hat C_1\otimes c_i) \oplus (\hat C_2 \otimes s_i))\big),$$
$$\hat y_i^{B}:\quad \hat A \oplus (\hat B \otimes f_i) \oplus (\hat C_1 \otimes g_i) \oplus (\hat C_2 \otimes h_i),$$

with $g_i = \mathrm{fl}(f_ic_i)$, $h_i = \mathrm{fl}(f_is_i)$. Algorithm A groups
as $f_i(C_1c_i)$; Algorithm B as $C_1(f_ic_i)$. Multiplication in floating point
is commutative but **not associative**, so the two need not agree.

Each is a backward-stable evaluation of the same expression in $O(1)$ operations,
so with $T_i = |\hat A| + |\hat Bf_i| + |\hat C_1g_i| + |\hat C_2h_i|$,

$$\big|\hat y_i^{A} - \hat y_i^{B}\big| \;\lesssim\; c\,u\,T_i$$

for a small constant $c$. Measured: $\max_i|\Delta\hat y_i| = 1.78\times10^{-15}$,
$\lVert\Delta\hat y\rVert_2 = 2.62\times10^{-14}$, and

$$\mathrm{rms}_i|\Delta\hat y_i| = 4.79\,u ,$$

against $\bar T = 8.82$ — i.e. $c \approx 0.5$ relative to $u\bar T$. The bound is
attained to within a small factor, as expected.

### Propagation to $S$

With $r = y - \hat y^{B}$ and $\Delta\hat y = \hat y^{A} - \hat y^{B}$,

$$S_A - S_B = \sum_i\Big[(r_i - \Delta\hat y_i)^2 - r_i^2\Big]
= -2\sum_i r_i\,\Delta\hat y_i + \lVert\Delta\hat y\rVert^2 .$$

The quadratic term is $\sim 10^{-27}$ and negligible, so the difference is
first-order in $\Delta\hat y$ and, crucially, **weighted by the residuals**:

| quantity | value |
|---|---|
| $S_A$ | $5.2867886392012755$ |
| $S_B$ | $5.2867886392012764$ |
| $S_A - S_B$ | $-8.88\times10^{-16}$ |
| relative | $1.68\times10^{-16}$ |
| first-order predictor $-2\sum_i r_i\Delta\hat y_i$ | $-1.11\times10^{-15}$ |

The predictor reproduces the observed difference to within 25%, confirming the
mechanism.

Two bounds bracket it. Cauchy–Schwarz gives the worst case,

$$\frac{|S_A - S_B|}{S} \le \frac{2\lVert r\rVert\,\lVert\Delta\hat y\rVert}{S}
= 2.28\times10^{-14},$$

which is $136\times$ above the observed value because it is tight only when
$\Delta\hat y \parallel r$. Treating the rounding errors as independent of $r$
gives the realistic estimate

$$\frac{|S_A - S_B|}{S} \approx \frac{2u\bar T}{\sqrt{S}} = 8.5\times10^{-16},$$

within a factor of five of the measured $1.68\times10^{-16}$.

**Conclusion.** The two implementations agree to about one unit in the last place
of $S$. The difference is not an approximation or a shortcut: it is the
unavoidable consequence of evaluating the same sum in a different association
order.

---

## 6. Conditioning — what the fusion does *not* change

For $X$ of full column rank with singular values $\sigma_1\ge\cdots\ge\sigma_4>0$,
$\kappa_2(X) = \sigma_1/\sigma_4$. From the SVD $X = U\Sigma V'$,

$$X'X = V\Sigma U'U\Sigma V' = V\Sigma^2V' \;\Longrightarrow\;
\sigma_i(X'X) = \sigma_i(X)^2 \;\Longrightarrow\; \kappa_2(X'X) = \kappa_2(X)^2 .$$

Verified: $\kappa_2(X) = 2.467\times10^{3}$ and $\kappa_2(X'X)$ computed directly
is $6.0875\times10^{6}$, matching $\kappa_2(X)^2$ to $1.5\times10^{-12}$ relative.

This is why `"qr"` exists: forming $X'X$ squares the condition number, costing
roughly $\log_{10}\kappa^2$ decimal digits, whereas a Householder QR of $X$ costs
about $\log_{10}\kappa$ when the residual is small. By Eckart–Young,
$1/\kappa_2(X)$ is the *relative* distance to the nearest rank-deficient matrix —
$4.05\times10^{-4}$ here, so the basis is a basis with four digits to spare.

Across the admissible box:

| $m$ | $\kappa_2(X)$ | $\kappa_2(X'X)$ | $\angle(\mathbf 1, \boldsymbol\tau^m)$ |
|---|---|---|---|
| 0.001 | $3.7\times10^{3}$ | $1.4\times10^{7}$ | 0.03° |
| 0.1 | $4.95\times10^{1}$ | $2.5\times10^{3}$ | 3.22° |
| 0.5 | $1.8\times10^{2}$ | $3.3\times10^{4}$ | 14.49° |
| 0.9 | $2.5\times10^{3}$ | $6.1\times10^{6}$ | 23.18° |

As $m\to0$, $\tau^m\to1$ and the second column rotates into the intercept. At the
enforced bound $m\ge0.1$ the normal equations cost about 3.4 of ~16 digits; even
at $m = 0.001$ the cost is ~7 digits. The regime where forming $X'X$ genuinely
fails needs $\kappa_2(X)\sim10^{8}$, i.e. $m\sim4\times10^{-8}$, far outside any
fit — but the dependency is real and is why `lower[1]` matters.

Since $\hat\beta$ is bitwise identical between Algorithms A and B, **none of this
is affected by the fusion.** It is a property of the `"crossprod"` method, before
and after.

---

## 7. Differences between `beta_method` values are a separate effect

Because $S$ is exactly quadratic in $\beta$ and $\nabla_\beta S(\hat\beta) = 0$,

$$S(\hat\beta + \delta) = S(\hat\beta) + \delta'(X'X)\delta$$

with **no first-order term** — an exact identity, not an expansion. A perturbation
of the coefficients is therefore damped quadratically in the objective, and the
damping depends on the direction of $\delta$ relative to the singular vectors of
$X$.

Measured against `"crossprod"`:

| method | $\lVert\delta\beta\rVert/\lVert\beta\rVert$ | observed $\Delta S/S$ | $\delta'(X'X)\delta / S$ |
|---|---|---|---|
| `"qr"` | $1.3\times10^{-15}$ | $3.0\times10^{-15}$ | $6.8\times10^{-27}$ |
| `"chol"` | $1.3\times10^{-15}$ | $3.5\times10^{-15}$ | $2.7\times10^{-26}$ |
| `"symengine"` | $2.4\times10^{-14}$ | $6.6\times10^{-15}$ | $2.5\times10^{-24}$ |

*This corrects an earlier characterisation.* For these four solvers the quadratic
term is $10^{-24}$ or smaller and explains none of the observed difference: the
$\Delta S$ between methods is dominated by rounding in the **evaluation** of $S$,
exactly as in §5, not by the difference in $\hat\beta$. The quadratic term becomes
the dominant contribution only for a solver whose $\hat\beta$ is much less
accurate — the hand-rolled Gauss–Jordan C++ prototype, at
$\lVert\delta\beta\rVert/\lVert\beta\rVert \approx 3.6\times10^{-10}$, is in that
regime.

---

## 8. Cost

Per objective evaluation at $n = 2428$, the fusion removes one full pass over the
series — one $\tau^m$, one $\log$, one $\cos$, one $\sin$ and the associated
temporaries:

| | µs per call |
|---|---|
| two-pass `"crossprod"` | 133 |
| fused `"crossprod"` | 90 |

Since the objective is evaluated 26,890 times per SPX F2 fit, this is the
difference between roughly 3.6 s and 2.4 s of objective evaluation.

---

## 9. Consequences for reproducibility

1. A calibration is bit-reproducible only against **the same `beta_method`**. The
   four methods differ at the $10^{-15}$ level in $S$ and at $10^{-15}$–$10^{-14}$
   in $\hat\beta$.
2. `"crossprod"` is additionally not bit-identical to a two-pass evaluation of the
   same normal equations, by $\approx 1$ ulp of $S$ (§5).
3. Neither difference moves a reported quantity. On the SPX fit, $t_c$ shifts in
   the 12th significant digit and $m$, $\omega$, the MPLE and all three likelihood
   intervals are unchanged.
4. Downstream figures are affected only through those trailing digits. Regenerating
   the thesis figure set with the fused objective left 35 of 38 PNGs byte-identical,
   the other three differing by 3–6 pixels at a maximum channel delta of 1/255.
5. `"symengine"` remains the path that reproduces the thesis calibration, and is
   deliberately left on the two-pass form so that it does.

---

## References

- G. H. Golub and V. Pereyra (1973), *The differentiation of pseudo-inverses and
  nonlinear least squares problems whose variables separate*, SIAM J. Numer. Anal.
  10(2), 413–432. — variable projection.
- N. J. Higham (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed.,
  SIAM. — chapters 19–20 for least-squares conditioning and the normal-equations
  vs QR comparison.
- G. H. Golub and C. F. Van Loan (2013), *Matrix Computations*, 4th ed., JHU Press.
  — SVD, Eckart–Young, perturbation of least squares.
- V. Filimonov and D. Sornette (2013), *A stable and robust calibration scheme of
  the log-periodic power law model*, Physica A 392(17), 3698–3707. — the
  separable parameterisation used here.
