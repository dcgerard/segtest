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
