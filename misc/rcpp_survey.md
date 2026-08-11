# Is there anything to be gained from Rcpp in lpplsF?

Survey and recommendation, 2026-08-10. Continues `misc/handoff_rcpp.txt`.

Everything below is measured on this machine (Apple Silicon, R 4.5.2, clang 17),
against the package as committed at `7d86146` ("Implement performance
enhancements"). The worktree's `R/` is byte-identical to dev's, so the
measurements are of dev's code; dev itself was never written to. Every timing is
the median of three repetitions, each in a **fresh R process**, on the real SPX
calibration (`n = 2428`, `fh = 500`, `num_searches = 1`, `factr = 1e7`).

---

## The short answer

**Yes, but far less than the component ratios implied — and the single largest
win in the package is not Rcpp at all.**

The `13.9x` and `33.4x` recorded for `compute_X_matrix()` and
`compute_H_matrix()` are not measuring "C++ beats R". They are measuring
**R's `for`-loop overhead against itself**. Both functions are trivially
vectorisable, and vectorising them in pure R recovers **7.5x and 10.4x** of
that — bit-identically, with no compiler.

Once the package is written the way R wants to be written, the whole remaining
Rcpp gain is about **1.2x**.

| SPX calibration | today | pure R, improved | Rcpp |
|---|---|---|---|
| F2 fit | 4.84 s | 3.62 s | 2.89 s |
| MPL fit | 8.55 s | 4.03 s | 3.39 s |
| MPL + the doc's plots | 10.11 s | 5.07 s | 4.19 s |
| rolling, K=21 MPL windows | 51.49 s | 29.72 s | 22.64 s |

Pure R gets **2.0x** on the doc's own SPX call. Rcpp gets **2.4x**.
The part only a compiler can buy is the last **1.21x**.

**Recommendation: do the pure-R work now; do not adopt `src/` yet.** Reasoning
in §6. A complete, tested, `R CMD check`-clean Rcpp build exists and is
preserved, so adopting later is a copy, not a project.

---

## 1. Method, and two traps worth recording

The survey timed each candidate on a realistic workload rather than in
isolation, and ran a control before attributing any difference to any change.
Two controls changed conclusions:

**(a) `pkgload::unload()` + `load_all()` in one session inflates the MPL stage
about 5x.** The first baseline measured the MPL stage at 19.1 s in a session
that had swapped package copies, and 3.7 s in one that had not — same code, same
data. Ruled out GC pressure explicitly (three scenarios: nothing retained, all
fits retained, 20,000 ballast objects live — 8.51 / 8.60 / 8.62 s, no effect).
The 19.1 s was an artifact of the reload. **Every number in this document comes
from a process that loaded exactly one package copy and then measured.**

**(b) `create_matrix_plot()` and `create_trace_plot()` take the objective as an
argument**, so a harness that builds its own closure measures the closure, not
the package. The first pass understated the compiled arm on both; §4 reports the
corrected numbers.

The `src/` arm was built so that **one binary serves both sides of every
comparison** — the compiled paths are selected at runtime by
`options(lpplsF.sse_engine=, lpplsF.mpl_engine=)`. The arm with `src/` present
but both engines set to `"r"` (`rr` below) reproduces the pure-R build to within
noise (F2 4.863 vs 4.887 s; MPL 8.667 vs 8.627 s), which is what licenses
attributing everything else to the engine rather than to the build.

The `Rprof()`-crashes-MPL bug from the handoff was not needed for this survey
(direct component timing worked), and has since been investigated separately —
see **`misc/rprof_mpl_crash.md`**. Short version: it is R's profiler colliding
with Apple Accelerate's `DGEMM`, reproducible in four lines of base R with no
packages loaded; not lpplsF, not symengine. `compute_mpl_loglik()` is the only
place in the package that multiplies two matrices, which is why F2 profiles fine
and MPL does not. Workaround for a profiling session:
`options(matprod = "internal")`.

---

## 2. Where the time actually goes

Objective-evaluation counts, from an instrumented build:

| workload | `SSE2` calls | `SSE3` | `beta_calculator` |
|---|---|---|---|
| SPX F2 | 26,890 | 11 | 27,402 |
| SPX MPL | 26,890 | 11 | 27,902 |
| SPX MPL + doc plots | 37,455 | 11 | 38,467 |
| synthetic (n=2000, fh=150) | 12,980 | 11 | 13,142 |

**Correction to the handoff.** The recorded "5,831 objective evaluations at
`factr = 1e7`" is `optim()`'s `counts["function"]`, which excludes the
finite-difference gradient evaluations. The true number of `SSE2` invocations is
**26,890** — 4.6x higher, and ~53.8 per candidate `tc`.

That correction overturns the handoff's conclusion that the F2 stage is "optim()
calls and not directly compilable". At 133 µs per call, those 26,890 evaluations
are **3.6 s of the 4.84 s F2 stage**. `optim()`'s own overhead is negligible;
L-BFGS-B is C already. **The F2 stage is almost entirely objective evaluation,
and the objective is a pure numeric function.**

---

## 3. The survey: every R-level loop over n or over a grid

| site | shape | runs when | cost (SPX) | verdict |
|---|---|---|---|---|
| `SSE2` objective, `lppls_fitter.R:376` | 26,890 calls/fit | **every fit** | 3.6 s | **fuse — the real target** |
| `compute_H_matrix`, `plotting.R:286` | `for` over n, 500/fit | every MPL fit | 0.56 s | **vectorise (10.4x, free)** |
| `compute_X_matrix`, `plotting.R:234` | `for` over n, 1000/fit | every MPL fit | 0.53 s | **vectorise (7.5x, free)** |
| `create_contour_data`, `plotting.R:446` | 101×101 objective grid | `cp`/`sp = TRUE` only | 1.37 s | follows the objective |
| `create_trace_plot` lattice, `plotting.R:668/709/760` | 100×100 per selector | `tp != 0` only | 1.38 s | follows the objective |
| `create_matrix_plot`, `plotting.R:600` | fh×3×3 = 4500 | `mp = TRUE` only | 0.68 s | follows the objective |
| trace-plot replay, `lppls_fitter.R:853` | O(k²) `optim()` calls | `tp != 0` only | **k = 15, 0.147 s** | **leave alone** |
| `lppls_rolling` window loop, `lppls_rolling.R:131` | K × `fit_lppls` | rolling | **no own overhead** | **nothing to compile** |
| `.sim_ar_garch`, `lr_trend_test.R:20,25` | 2 serial recursions | B× in bootstrap | **0.13 s per 199 paths** | **leave alone** |
| `eval_lppls` mode 0 | vectorised already | everywhere | 51 µs | subsumed by the fusion |
| per-`tc` tibble/dplyr block | 500 × arrange/filter/rbind | every F2 fit | **0.62 s** | R-level, not Rcpp |

Three of the handoff's open items resolve as **non-candidates**:

- **The O(k²) trace-plot replay is not a hot spot.** `k` is the L-BFGS-B function
  count for the best search — **15** on SPX — so the replay is 15 `optim()` calls
  totalling 0.147 s. The quadratic shape is real but the constant never gets
  large at usable settings.
- **`lppls_rolling()` has no overhead of its own.** K=5 windows through
  `lppls_rolling()` took 8.20 s; the same five fits called directly took 8.79 s.
  The loop is a driver; all its time is `fit_lppls`.
- **`.sim_ar_garch()` is 10x in C++ and it does not matter.** 199 paths cost
  0.13 s, inside a bootstrap whose K×B rolling calibrations cost ~7,300 s. It is
  0.002% of the run. `lr_trend_test()` gets faster only because the *fits* inside
  it get faster.

One item is a **pure-R opportunity Rcpp does not touch**: the per-`tc`
bookkeeping block (`tibble` row assignment, `dplyr::arrange`/`filter`, growing
`rbind`) costs **0.62 s per fit** — 1,242 µs per candidate `tc`, measured
identically in both arms. That is 13% of today's F2 stage and 21% of the
compiled one. Preallocated numeric vectors would remove most of it.

---

## 4. Per-call numbers (n = 2428, µs/call)

| routine | today | vectorised / fused **R** | **C++** | R gain | extra from C++ |
|---|---|---|---|---|---|
| `SSE2` objective | 133.0 | 90.0 | **36.5** | 1.5x | 2.5x |
| `compute_X_matrix` | 526.0 | **70.0** | 42.0 | **7.5x** | 1.7x |
| `compute_H_matrix` | 1122.0 | **108.0** | 34.0 | **10.4x** | 3.2x |
| beta solver | 77.5 | — | 33.0 | — | 2.3x |
| `eval_lppls` (mode 0) | 51.0 | — | 33.5 | — | 1.5x |
| `.sim_ar_garch` (n=2000) | 600.0 | — | 60.0 | — | 10x |

Grid scans, with the objective correctly passed in:

| routine | today | C++ |
|---|---|---|
| `create_contour_data()` 101×101 | 1.374 s | 0.609 s |
| `create_trace_plot()` selector 1 | 1.381 s | 0.695 s |
| `create_matrix_plot()` 4500 cells | 0.683 s | 0.395 s |

**This table is the whole argument.** For the two MPL matrices, C++ is buying
1.7x and 3.2x over vectorised R — not 12x and 33x. The large ratios in the
handoff compared C++ against an *interpreted `for` loop*, which is a fixable
property of the R code rather than a property of R.

For the objective it is the other way round: pure R can only get 1.5x (by not
building the `tau^m`/`cos`/`sin` basis twice — once in `design()`, once in
`eval_lppls()`), and C++ gets 3.6x. **The objective is the one place where a
compiler earns something R cannot.**

---

## 5. Whole-fit results and invariance

Median of 3, fresh process each. `rr` is the control described in §1.

| arm | F2 | MPL | MPL + doc plots |
|---|---|---|---|
| today, pure R | 4.838 | 8.545 | 10.106 |
| `src/` present, engines = R (**control**) | 4.863 | 8.667 | — |
| vectorised MPL loops only | 4.775 | 5.252 | 6.790 |
| + fused pure-R objective | 3.615 | 4.028 | 5.071 |
| compiled objective only | 2.934 | 6.618 | — |
| compiled MPL matrices only | 4.861 | 5.381 | — |
| **Rcpp, both** | **2.894** | **3.394** | **4.191** |

Rolling, K = 21 MPL windows (n = 2000, fh = 150):
today **51.49 s** (2.45 s/window) → pure R improved **29.72 s** (1.42) → Rcpp
**22.64 s** (1.08).

### Invariance

- **Vectorising `compute_X_matrix`/`compute_H_matrix` is bit-identical.** Max
  absolute difference 0.000e+00 against the loops; the full SPX F2, MPL and
  doc-plot digests (tc, m, ω, A, B, C1, C2, D, value, MPLE, MLL, all three
  likelihood intervals, 455 defined points) are identical to the current
  baseline. **No figure needs re-rendering.**
- **Changing the objective is not bit-identical — in R or in C++.** Both fused
  forms shift `tc` in the 12th significant digit
  (2717.00880440869 → …40381 in R, → …42704 in C++; ≤ 7e-12 relative). In every
  case `m`, `ω`, the MPLE (2717), all three likelihood intervals
  ([2696,2757]/[2696,2757]/[2696,2744]) and the 455 defined points are unchanged.
  Over a 21-window rolling calibration, `tc` moves by ≤ 2.0e-7 absolute
  (1.0e-10 relative) and **all six likelihood-interval columns and `tc_hat_mpl`
  are bit-identical**. This is the same character of change as the `factr` fix,
  and it is a property of fusing the objective, not of C++.
- **Tests: 556 pass / 0 fail** on the pure-R improvement, and on the compiled
  build under both engines — matching the current package exactly.

---

## 6. The cost of `src/`, measured

| | pure R | with `src/` |
|---|---|---|
| `R CMD check --as-cran` | 2 NOTEs | **2 NOTEs, same ones** |
| `checking compiled code` | — | **OK** |
| source tarball | 163,658 B | 170,142 B (+6.5 KB) |
| `load_all()` warm | 0.62–0.65 s | 0.78 s |
| `load_all()` cold (after editing `.cpp`) | — | 2.90 s |
| test suite | 556 / 0 | 556 / 0 |

Adding `src/` introduced **no new check finding**. (An earlier run showed a
WARNING and an extra NOTE; both were my own patching artifacts — an
un-regenerated `.Rd` after changing `beta_method`'s formals, and a `.Rbuildignore`
that `list.files()` had not copied. With both fixed the two checks are identical.)

The handoff's point stands: **`symengine` is a compiled package**, so a source
install of lpplsF already needs a C++ toolchain, and `symengine` is by far the
heavier compile. The toolchain argument no longer counts against `src/`. For
CRAN binary users it was never a factor either way.

### So why not adopt

1. **The measured incremental gain is 1.19–1.31x**, once the pure-R work is
   done. That is a real number but a small one to buy with a permanent
   structural change.
2. **Most of the apparent gain is free.** Vectorising two loops takes the MPL fit
   from 8.55 s to 5.25 s, bit-identically, today, with no compiler and no
   re-render. That should happen regardless of what is decided about `src/`.
3. **It adds a correctness surface.** The prototype's hand-rolled Gauss-Jordan
   is *less* accurate than R's `solve()` on the same normal equations
   (3.6e-10 vs the `crossprod` path). Harmless inside `0.1 ≤ m ≤ 0.9`, but it is
   hand-written linear algebra to maintain, and using RcppArmadillo instead would
   trade that for a much heavier dependency.
4. **Sequencing.** The stated plan is thesis first, packaging decisions last, and
   this is a first CRAN submission by a solo maintainer — `src/` multiplies the
   platforms that have to build cleanly at exactly the wrong moment.
5. **The better axis is still the sweep.** Parallelising the outer sweep gives
   far more than 1.2x and needs no toolchain, consistent with the earlier
   decision to decline parallelism *inside* `fit_lppls()`.

### What would change the answer

If the sweep layer becomes the dominant cost — thousands of rolling calibrations
— then 2.27x versus 1.73x is the difference between roughly 10 hours and 7.6
hours per sweep, and `src/` starts paying for itself. That is a roadmap
judgement, not a measurement one. Everything needed is built and verified.

---

## 7. If `src/` is adopted anyway: the API

Following the `beta_method` pattern — implementations shipped as options, best
performing as default, never a silent replacement. Two knobs, because the two
compiled pieces have **different invariance properties** and should not be
forced to move together:

```r
fit_lppls(...,
          beta_method = c("cpp", "crossprod", "qr", "chol", "symengine"),
          mpl_engine  = c("cpp", "r"))
```

- `beta_method = "cpp"` selects the fused compiled objective *and* solver. It
  belongs on this argument because it is genuinely a fifth answer to "how are the
  linear parameters solved" — but the docs must say that it also replaces the
  SSE objective, since that is both what makes it fast and why it is not
  bit-identical to `"crossprod"`.
- `mpl_engine` is separate because the compiled MPL matrices **are** bit-identical,
  so that switch is safe unconditionally, whereas the objective switch is a
  12th-digit change a user may not want.
- `"symengine"` stays as the thesis-checkable path, as decided.

The working build lives at
`$SCRATCH/pkg_rcpp` (see §8); `misc/rcpp_survey_prototypes.cpp` holds the
sources, including a fused `mpl_moments_cpp()` that accumulates `X'X`,
`X_hat'X_tc` and `H` in one pass without materialising either n×6 matrix, and
leaves `det()` to R so the `det1 <= 0` guard stays on LAPACK.

---

## 8. Recommended order of work

1. **Vectorise `compute_X_matrix()` and `compute_H_matrix()`** (`R/plotting.R`).
   Free, bit-identical, 1.63x on an MPL fit. Code in §9. *No cache regeneration
   or figure re-render needed — output is unchanged.*
2. **Flatten the per-`tc` bookkeeping** in the F2 loop — 0.62 s/fit, pure R.
3. **Optionally fuse the objective in R** (1.34x on F2). Needs the same
   verification the `factr` change got: recompute the three single fits, confirm
   reported quantities, control-render the figures.
4. **Revisit `src/` when the sweep layer exists.**

Steps 1–3 give the 2.0x column. Step 4 is the remaining 1.2x.

## 9. The Tier-1 patch (verified bit-identical, 556 tests pass)

`dev` is read-only, so this is here rather than applied. Replaces the bodies of
`compute_X_matrix()` and `compute_H_matrix()` in `R/plotting.R`; signatures,
roxygen and call sites are unchanged.

```r
compute_X_matrix <- function(Psi, tc, t) {
  m <- Psi[1]; omega <- Psi[2]; b <- Psi[4]; c1 <- Psi[5]; c2 <- Psi[6]
  tau <- tc - t; log_tau <- log(tau); tau_m <- tau^m
  cos_term <- cos(omega * log_tau); sin_term <- sin(omega * log_tau)
  cbind(tau_m * log_tau * (b + c1 * cos_term + c2 * sin_term),
        tau_m * log_tau * (-c1 * sin_term + c2 * cos_term),
        1, tau_m, tau_m * cos_term, tau_m * sin_term)
}

compute_H_matrix <- function(Psi, tc, log_p, t) {
  m <- Psi[1]; omega <- Psi[2]; a <- Psi[3]
  b <- Psi[4]; c1 <- Psi[5]; c2 <- Psi[6]
  tau <- tc - t; log_tau <- log(tau); tau_m <- tau^m
  cos_term <- cos(omega * log_tau); sin_term <- sin(omega * log_tau)
  lt2 <- log_tau * log_tau
  res <- log_p - (a + tau_m * (b + c1 * cos_term + c2 * sin_term))
  H <- matrix(0, nrow = 6, ncol = 6)
  H[1, 1] <- sum(res * tau_m * lt2 * (b + c1 * cos_term + c2 * sin_term))
  H[1, 2] <- sum(res * tau_m * lt2 * (-c1 * sin_term + c2 * cos_term))
  H[1, 4] <- sum(res * tau_m * log_tau)
  H[1, 5] <- sum(res * tau_m * log_tau * cos_term)
  H[1, 6] <- sum(res * tau_m * log_tau * sin_term)
  H[2, 2] <- sum(res * (-1) * tau_m * lt2 * (c1 * cos_term + c2 * sin_term))
  H[2, 5] <- sum(res * (-1) * tau_m * log_tau * sin_term)
  H[2, 6] <- sum(res * tau_m * log_tau * cos_term)
  ## Symmetry -- H holds residual-weighted second derivatives, so both
  ## triangles must be filled or det(X'X - H) is not the observed information.
  H[2, 1] <- H[1, 2]; H[4, 1] <- H[1, 4]; H[5, 1] <- H[1, 5]
  H[6, 1] <- H[1, 6]; H[5, 2] <- H[2, 5]; H[6, 2] <- H[2, 6]
  H
}
```

Note the vectorised `compute_X_matrix()` relies on `cbind()` recycling the
scalar `1` into the third column, which is what the loop's `X[i, 3] <- 1` did.

---

## 10. Artifacts

In this worktree:
- `misc/rcpp_survey.md` — this document
- `misc/rcpp_survey_prototypes.cpp` — all prototypes, verified (agreement in §4;
  the rejected one-pass `y'y - b'X'y` objective is kept and labelled, since its
  8.7e-11 error is the reason two passes are used)

In the session scratchpad (`.../scratchpad/`), not in any repo:
`pkg_control` (pure-R copy of dev), `pkg_rvec` (pure-R improvement),
`pkg_rcpp` (the `src/` build, `R CMD check` clean), `pkg_instr` (call counters),
and scripts `01`–`12` with their logs.
