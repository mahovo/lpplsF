# Tests for residual_qq_plot(): residual QQ diagnostic against a normal or a
# fitted Student-t reference distribution, with green reference lines at the
# chosen quantiles.

test_that("residual_qq_plot (normal) is a ggplot with qq points and one h/v ref line each", {
  set.seed(1)
  p <- residual_qq_plot(rnorm(200))
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomPoint" %in% geoms)          # stat_qq points
  expect_equal(sum(geoms == "GeomHline"), 1L)  # sample-quantile refs
  expect_equal(sum(geoms == "GeomVline"), 1L)  # theoretical-quantile refs
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("normal reference lines sit at the sample and normal-theory quantiles", {
  set.seed(2)
  r  <- rnorm(500)
  p  <- residual_qq_plot(r, ref = c(0.1, 0.9))
  hl <- Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)[[1]]
  vl <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)[[1]]
  expect_equal(sort(unname(hl$data$yintercept)),
               sort(unname(stats::quantile(r, c(0.1, 0.9)))))
  expect_equal(sort(unname(vl$data$xintercept)), sort(stats::qnorm(c(0.1, 0.9))))
})

test_that("distribution = 't' builds silently", {
  skip_if_not_installed("MASS")
  set.seed(3)
  expect_silent(ggplot2::ggplot_build(residual_qq_plot(rt(400, df = 5),
                                                       distribution = "t")))
})

test_that("residual_qq_plot rejects non-numeric or too-short input", {
  expect_error(residual_qq_plot("a"))
  expect_error(residual_qq_plot(1))
})
