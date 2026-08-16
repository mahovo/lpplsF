## The four linear-parameter solvers, and a computer-algebra oracle for the
## hand-transcribed MPL derivative matrices.

METHODS <- c("crossprod", "qr", "chol", "symengine")

sample_case <- function(n = 400, m = 0.6, omega = 8, noise = 0.01, seed = 1) {
  set.seed(seed)
  t  <- seq_len(n)
  tc <- n + 21
  lp <- eval_lppls(t, 8, -0.025, 0.0015, 0.001, tc, m, omega) +
    stats::rnorm(n) * noise
  list(t = t, log_p = lp, tc = tc, m = m, omega = omega)
}


test_that("every solver returns the same beta, shaped the same way", {
  d <- sample_case()
  out <- lapply(METHODS, function(mm) {
    create_beta_calculator(mm)(d$log_p, d$t, d$tc, d$m, d$omega)
  })
  names(out) <- METHODS

  for (mm in METHODS) {
    expect_length(out[[mm]], 4)
    expect_named(out[[mm]], c("A", "B", "C1", "C2"))
    expect_true(all(is.finite(out[[mm]])))
  }
  ## The solvers agree far more tightly than anything the model interprets.
  for (mm in setdiff(METHODS, "crossprod")) {
    expect_equal(unname(out[[mm]]), unname(out[["crossprod"]]), tolerance = 1e-8)
  }
})


test_that("the solvers recover a known beta exactly", {
  ## y = X beta with no noise, so beta is known rather than merely agreed on.
  d <- sample_case(noise = 0)
  beta_true <- c(A = 8, B = -0.025, C1 = 0.0015, C2 = 0.001)
  for (mm in METHODS) {
    got <- create_beta_calculator(mm)(d$log_p, d$t, d$tc, d$m, d$omega)
    expect_equal(unname(got), unname(beta_true), tolerance = 1e-6,
                 info = paste("method:", mm))
  }
})


test_that("create_beta_calculator validates its method", {
  expect_error(create_beta_calculator("nonsense"))
  ## Default is the documented one.
  expect_equal(create_beta_calculator()(1:10 + 0, 1:10, 31, 0.6, 8),
               create_beta_calculator("crossprod")(1:10 + 0, 1:10, 31, 0.6, 8))
})


test_that("fit_lppls accepts beta_method and records it", {
  d <- sample_case(n = 200)
  f <- fit_lppls(log_price = d$log_p, fh = 10, hold_out = 0, tc_init = 205,
                 m_init = 0.5, o_init = 8, num_searches = 1, mode = "F2",
                 beta_method = "qr")
  expect_identical(f$fit_args$beta_method, "qr")
  expect_error(
    fit_lppls(log_price = d$log_p, fh = 10, hold_out = 0, tc_init = 205,
              mode = "F2", beta_method = "nonsense"))
})


test_that("the choice of solver does not move the estimate", {
  d <- sample_case(n = 200)
  fits <- lapply(METHODS, function(mm) {
    fit_lppls(log_price = d$log_p, fh = 10, hold_out = 0, tc_init = 205,
              m_init = 0.5, o_init = 8, num_searches = 1, mode = "F2",
              beta_method = mm)$fit[[1]]
  })
  names(fits) <- METHODS
  for (mm in setdiff(METHODS, "crossprod")) {
    for (k in c("tc", "m", "omega", "value")) {
      expect_equal(fits[[mm]][[k]], fits[["crossprod"]][[k]], tolerance = 1e-4,
                   info = paste(mm, k))
    }
  }
})


test_that("MPL mode works under every solver and agrees on the intervals", {
  d <- sample_case(n = 200)
  run <- function(mm) {
    suppressWarnings(
      fit_lppls(log_price = d$log_p, fh = 20, hold_out = 0, tc_init = 210,
                m_init = 0.5, o_init = 8, num_searches = 1, mode = "MPL",
                mpl_cutoff = c(0.05, 0.5), beta_method = mm))
  }
  ref <- run("crossprod")$mpl_output
  for (mm in setdiff(METHODS, "crossprod")) {
    got <- run(mm)$mpl_output
    expect_equal(got$tc_hat_mpl, ref$tc_hat_mpl, info = mm)
    expect_equal(got$LI, ref$LI, info = mm)
  }
})


## ---- computer-algebra oracle -------------------------------------------------
## compute_X_matrix() and compute_H_matrix() transcribe Filimonov (2017)
## eq. (36) and (37) by hand. Deriving them symbolically from the model gives an
## independent check: a transcription slip in either shows up here. The
## asymmetric-H defect was exactly such a slip.

lppls_symbolic <- function() {
  tau <- symengine::Symbol("tau"); m <- symengine::Symbol("m")
  w   <- symengine::Symbol("w");   A <- symengine::Symbol("A")
  B   <- symengine::Symbol("B");   C1 <- symengine::Symbol("C1")
  C2  <- symengine::Symbol("C2")
  f <- A + tau^m * (B + C1 * cos(w * log(tau)) + C2 * sin(w * log(tau)))
  list(f = f, args = c(tau, m, w, A, B, C1, C2),
       pars = list(m = m, w = w, A = A, B = B, C1 = C1, C2 = C2))
}

## Evaluate a symbolic expression over the sample, at fixed parameters.
sym_eval <- function(expr, sy, taus, Psi) {
  vis <- symengine::DoubleVisitor(expr, args = sy$args)
  n <- length(taus)
  vis(taus, rep(Psi[1], n), rep(Psi[2], n), rep(Psi[3], n),
      rep(Psi[4], n), rep(Psi[5], n), rep(Psi[6], n))
}


test_that("compute_X_matrix matches symbolic differentiation of the model", {
  skip_if_not_installed("symengine")
  sy <- lppls_symbolic()
  n <- 150; t <- seq_len(n); tc <- n + 21
  Psi <- c(0.6, 8, 8, -0.025, 0.0015, 0.001)   # (m, omega, A, B, C1, C2)
  taus <- tc - t

  ## Column order in the hand-written matrix is (m, omega, A, B, C1, C2).
  cols <- c("m", "w", "A", "B", "C1", "C2")
  X_cas <- vapply(cols, function(nm) {
    sym_eval(symengine::D(sy$f, sy$pars[[nm]]), sy, taus, Psi)
  }, numeric(n))

  expect_equal(unname(compute_X_matrix(Psi, tc, t)), unname(X_cas),
               tolerance = 1e-10)
})


test_that("compute_H_matrix matches symbolic second derivatives and is symmetric", {
  skip_if_not_installed("symengine")
  sy <- lppls_symbolic()
  n <- 150; t <- seq_len(n); tc <- n + 21
  Psi <- c(0.6, 8, 8, -0.025, 0.0015, 0.001)
  set.seed(2)
  log_p <- eval_lppls(t, Psi[3], Psi[4], Psi[5], Psi[6], tc, Psi[1], Psi[2]) +
    stats::rnorm(n) * 0.01
  taus <- tc - t

  ## H[j, k] = sum_i residual_i * d2 f_i / d psi_j d psi_k.
  fitted <- eval_lppls(t, Psi[3], Psi[4], Psi[5], Psi[6], tc, Psi[1], Psi[2])
  res <- log_p - fitted
  cols <- c("m", "w", "A", "B", "C1", "C2")
  H_cas <- matrix(0, 6, 6)
  for (j in seq_along(cols)) {
    for (k in seq_along(cols)) {
      d2 <- symengine::D(symengine::D(sy$f, sy$pars[[cols[j]]]), sy$pars[[cols[k]]])
      if (as.character(d2) == "0") next
      H_cas[j, k] <- sum(res * sym_eval(d2, sy, taus, Psi))
    }
  }

  H_hand <- compute_H_matrix(Psi, tc, log_p, t)
  expect_equal(unname(H_hand), unname(H_cas), tolerance = 1e-8)

  ## Mixed partials commute, so H must be symmetric. This is the check that the
  ## hand-written version previously failed.
  expect_equal(H_hand, t(H_hand))
  expect_true(isSymmetric(unname(H_hand)))
})


test_that("the derivative entries the hand-written H omits are genuinely zero", {
  skip_if_not_installed("symengine")
  sy <- lppls_symbolic()
  ## Everything among the linear parameters, and A against anything, must vanish:
  ## the model is linear in A, B, C1, C2.
  for (pair in list(c("A", "A"), c("A", "m"), c("A", "w"), c("B", "B"),
                    c("B", "w"), c("B", "C1"), c("C1", "C1"), c("C1", "C2"))) {
    d2 <- symengine::D(symengine::D(sy$f, sy$pars[[pair[1]]]), sy$pars[[pair[2]]])
    expect_equal(as.character(d2), "0",
                 info = paste("d2 /", pair[1], pair[2]))
  }
})


# ---------------------------------------------------------------------------
# The objective. create_sse_calculator() fuses the "crossprod" path -- the
# design matrix is built once and the residuals come from it -- while every
# other solver keeps the original two-pass evaluation.
# ---------------------------------------------------------------------------

sse_two_pass <- function(method) {
  bc <- create_beta_calculator(method)
  function(log_p, t, tc, m, omega) {
    b <- bc(log_p, t, tc, m, omega)
    sum((log_p - eval_lppls(t, b["A"], b["B"], b["C1"], b["C2"],
                            tc, m, omega, mode = 0))^2, na.rm = TRUE)
  }
}

test_that("symengine stays on the two-pass form, bit-identically", {
  d <- sample_case()
  ## This is the guarantee that lets "symengine" reproduce the thesis
  ## calibration exactly rather than merely to rounding.
  expect_identical(
    lpplsF:::create_sse_calculator("symengine")(d$log_p, d$t, d$tc, d$m, d$omega),
    sse_two_pass("symengine")(d$log_p, d$t, d$tc, d$m, d$omega))
})

test_that("the fused objectives match their two-pass forms to machine precision", {
  d <- sample_case()
  for (mth in c("crossprod", "qr", "chol")) {
    a <- lpplsF:::create_sse_calculator(mth)(d$log_p, d$t, d$tc, d$m, d$omega)
    b <- sse_two_pass(mth)(d$log_p, d$t, d$tc, d$m, d$omega)
    ## Same quantity summed in a different order: equal to within rounding,
    ## and deliberately not asserted bit-identical.
    expect_equal(a, b, tolerance = 1e-12, info = mth)
  }
})

test_that("all four objectives agree across the admissible box", {
  d <- sample_case()
  grid <- expand.grid(m = c(0.1, 0.35, 0.6, 0.9), omega = c(6, 9.5, 13))
  ref <- lpplsF:::create_sse_calculator("crossprod")
  for (mth in c("qr", "chol", "symengine")) {
    f <- lpplsF:::create_sse_calculator(mth)
    for (i in seq_len(nrow(grid))) {
      expect_equal(f(d$log_p, d$t, d$tc, grid$m[i], grid$omega[i]),
                   ref(d$log_p, d$t, d$tc, grid$m[i], grid$omega[i]),
                   tolerance = 1e-8,
                   info = sprintf("%s at m=%g omega=%g", mth, grid$m[i], grid$omega[i]))
    }
  }
})

test_that("the objective equals the residual sum of squares at the returned beta", {
  d <- sample_case()
  for (mth in c("crossprod", "qr", "chol", "symengine")) {
    b <- create_beta_calculator(mth)(d$log_p, d$t, d$tc, d$m, d$omega)
    fitted <- eval_lppls(d$t, b["A"], b["B"], b["C1"], b["C2"],
                         d$tc, d$m, d$omega, mode = 0)
    expect_equal(lpplsF:::create_sse_calculator(mth)(d$log_p, d$t, d$tc, d$m, d$omega),
                 sum((d$log_p - fitted)^2), tolerance = 1e-10, info = mth)
  }
})

test_that("create_sse_calculator validates its method", {
  expect_error(lpplsF:::create_sse_calculator("nope"), "should be one of")
})
