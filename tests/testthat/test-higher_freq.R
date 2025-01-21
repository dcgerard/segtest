test_that("expit and logit are inverses", {
  expect_equal(logit(expit(c(0.3, 0.4))), c(0.3, 0.4))
})

testthat("real_to_simplex and simplex_to_real are inverses", {
  expect_equal(simplex_to_real(real_to_simplex(c(1, 3))), c(1, 3))
})

testthat("offspring_gf_2 and f1_gf_dr are same for tetraploids", {
  expect_equal(
    offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 1),
    f1_gf_dr(alpha = 1/12, g1 = 2, g2 = 1, ploidy = 4)
  )

  expect_equal(
    offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 1, p2 = 3),
    f1_gf_dr(alpha = 1/12, g1 = 1, g2 = 3, ploidy = 4)
  )
})
