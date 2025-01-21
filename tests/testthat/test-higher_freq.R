test_that("expit and logit are inverses", {
  expect_equal(logit(expit(c(0.3, 0.4))), c(0.3, 0.4))
})

test_that("real_to_simplex and simplex_to_real are inverses", {
  expect_equal(simplex_to_real(real_to_simplex(c(1, 3))), c(1, 3))
})

test_that("offspring_gf_2 and f1_gf_dr are same for tetraploids", {
  expect_equal(
    offspring_gf_2(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 1),
    f1_gf_dr(alpha = 1/12, g1 = 2, g2 = 1, ploidy = 4)
  )

  expect_equal(
    pvec_tet_2(alpha = 1/12, xi = 1/3, ell = 0),
    gamfreq_dr(alpha = 1/12, g = 0, ploidy = 4, log_p = FALSE)
  )

  expect_equal(
    pvec_tet_2(alpha = 1/12, xi = 1/3, ell = 1),
    gamfreq_dr(alpha = 1/12, g = 1, ploidy = 4, log_p = FALSE)
  )

  expect_equal(
    pvec_tet_2(alpha = 1/12, xi = 1/3, ell = 2),
    gamfreq_dr(alpha = 1/12, g = 2, ploidy = 4, log_p = FALSE)
  )

  expect_equal(
    pvec_tet_2(alpha = 1/12, xi = 1/3, ell = 3),
    gamfreq_dr(alpha = 1/12, g = 3, ploidy = 4, log_p = FALSE)
  )

  expect_equal(
    pvec_tet_2(alpha = 1/12, xi = 1/3, ell = 4),
    gamfreq_dr(alpha = 1/12, g = 4, ploidy = 4, log_p = FALSE)
  )
})
