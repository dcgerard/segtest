test_that("par_to_gf works", {
  rule <- list(
    list(ploidy = 4, g = 0, type = "mix"),
    list(ploidy = 8, g = 8, type = "polysomic"),
    list(outlier = FALSE)
    )
  par <- c()
  expect_equal(par_to_gf(par = par, rule = rule), c(0, 0, 0, 0, 1, 0, 0), tolerance = 1e-3)

  rule <- list(
    list(ploidy = 4, g = 2, type = "mix"),
    list(ploidy = 8, g = 8, type = "polysomic"),
    list(outlier = FALSE)
    )
  par <- c(-5)
  expect_no_error(par_to_gf(par = par, rule = rule))

})

test_that("par_to_gam and gam_to_par are inverses", {
  rule <- list(
    list(ploidy = 4, g = 2, type = "mix"),
    list(ploidy = 8, g = 6, type = "polysomic"),
    list(outlier = TRUE)
    )
  par <- c(-2, 0.1, 0.05, 0.03)

  ret <- gam_to_par(par_to_gam(par = par, rule = rule))
  expect_equal(par, ret$par)
  expect_equal(rule, ret$rule)

  rule <- list(
    list(ploidy = 4, g = 1, type = "mix_dr"),
    list(ploidy = 8, g = 0, type = "polysomic"),
    list(outlier = TRUE)
    )
  par <- c(0.1, 0.2)

  ret <- gam_to_par(par_to_gam(par = par, rule = rule))
  expect_equal(par, ret$par)
  expect_equal(rule, ret$rule)

  rule <- list(
    list(ploidy = 12, g = 6, type = "mix"),
    list(ploidy = 12, g = 6, type = "polysomic"),
    list(outlier = TRUE)
    )
  par <- c(-2, -1, 2, 0.1, 0.05, 0.01, 0.2)

  ret <- gam_to_par(par_to_gam(par = par, rule = rule))
  expect_equal(par, ret$par)
  expect_equal(rule, ret$rule)
})


test_that("fixed parameterizations work", {
  lower_val <- -20
  upper_val <- 20
  TOL <- sqrt(.Machine$double.eps)

  ## fix gamma
  rule <- list(
    list(ploidy = 12, g = 6, type = "mix", gamma = c(0.1, 0.5, 0.2, 0.2)),
    list(ploidy = 12, g = 6, type = "polysomic"),
    list(outlier = TRUE)
    )
  par <- c(0.1, 0.05, 0.01, 0.2)
  fix_list <- list(
    list(alpha = FALSE, gamma = TRUE, beta = FALSE),
    list(alpha = FALSE, gamma = FALSE, beta = FALSE),
    list(pi = FALSE)
  )
  expect_equal(par_to_gam(par = par, rule = rule)[[1]]$gamma, c(0.1, 0.5, 0.2, 0.2))
  ret <- gam_to_par(gam = par_to_gam(par = par, rule = rule), fix_list = fix_list, ob = 0.03)
  expect_equal(ret$par, par)
  expect_equal(ret$rule, rule)
  expect_equal(ret$lower, rep(TOL, 4))
  expect_equal(ret$upper, c(drbounds(ploidy = 12, model = "ces"), 0.03))

  ## fix alpha
  rule <- list(
    list(ploidy = 12, g = 6, type = "mix"),
    list(ploidy = 12, g = 6, type = "polysomic", alpha = c(0.2, 0.1, 0.05)),
    list(outlier = TRUE)
    )
  par <- c(-2, -1, 2, 0.2)
  fix_list <- list(
    list(alpha = FALSE, gamma = FALSE, beta = FALSE),
    list(alpha = TRUE, gamma = FALSE, beta = FALSE),
    list(pi = FALSE)
  )
  expect_equal(par_to_gam(par = par, rule = rule)[[2]]$alpha, c(0.2, 0.1, 0.05))

  ret <- gam_to_par(gam = par_to_gam(par = par, rule = rule), fix_list = fix_list, ob = 0.03)
  expect_equal(ret$par, par)
  expect_equal(ret$rule, rule)
  expect_equal(ret$lower, c(rep(lower_val, times = 3), TOL))
  expect_equal(ret$upper, c(rep(upper_val, times = 3), 0.03))

  ## fix outlier
  rule <- list(
    list(ploidy = 12, g = 6, type = "mix"),
    list(ploidy = 12, g = 6, type = "polysomic"),
    list(outlier = TRUE, pi = 0.5)
    )
  par <- c(1, 2, 1, -2, -1, 2)
  fix_list <- list(
    list(alpha = FALSE, gamma = FALSE, beta = FALSE),
    list(alpha = FALSE, gamma = FALSE, beta = FALSE),
    list(pi = TRUE)
  )
  expect_equal(par_to_gam(par = par, rule = rule)[[3]]$pi, 0.5)

  ret <- gam_to_par(gam = par_to_gam(par = par, rule = rule), fix_list = fix_list, ob = 0.03)
  expect_equal(ret$par, par)
  expect_equal(ret$rule, rule)
  expect_equal(ret$lower, c(rep(lower_val, 3), rep(TOL, 3)))
  expect_equal(ret$upper, c(rep(upper_val, 3), drbounds(ploidy = 12, model = "ces")))

})

