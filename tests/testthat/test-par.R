test_that("par_to_gf works", {
  rule <- list(
    list(ploidy = 4, g = 0, type = "mix"),
    list(ploidy = 8, g = 8, type = "polysomic"),
    list(outlier = FALSE)
    )
  par <- c()
  expect_equal(par_to_gf(par = par, rule = rule), c(0, 0, 0, 0, 1, 0, 0))

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
