test_that("hard cases work", {
  x <- c(6L, 87L, 429L, 392L, 86L)
  g1 <- 3
  g2 <- 2
  par = c(1e-5, 1e-5, 1 - 1e-5)
  expect_no_error(
    obj_dr_pp(par = par, x = x, g1 = g1, g2 = g2)
  )
})
