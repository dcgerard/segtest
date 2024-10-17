test_that("impossible genotypes calculated correctly", {
  xmat <- structure(
    c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
      4, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
      3, 4, 4, 5, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 0, 0,
      0, 0, 1, 1, 1, 2, 2, 3, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 0, 0, 0, 0, 1, 1, 1, 2, 2,
      3, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3,
      0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0,
      0, 1, 0, 0, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1,
      2, 0, 1, 0, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 0, 1,
      2, 3, 0, 1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 1, 2,
      3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 0, 1, 2, 3, 0, 1, 2, 0, 1,
      0, 0, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0,
      0, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0,
      1, 0, 0, 0, 5, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1,
      0, 1, 0, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0, 3, 2,
      1, 0, 2, 1, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 4, 3, 2,
      1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0,
      0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0,
      2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0), dim = c(126L, 5L), dimnames = list(NULL, NULL))

  pdf <- expand.grid(p1 = 0:4, p2 = 0:4, dr = c(TRUE, FALSE))

  TOL <- sqrt(.Machine$double.eps)
  for (i in seq_len(nrow(pdf))) {
    gf <- offspring_gf_2(
      alpha = ifelse(pdf$dr[[i]], 1/6, 0),
      xi1 = 1/3,
      xi2 = 1/3,
      p1 = pdf$p1[[i]],
      p2 = pdf$p2[[i]])

    if (any(gf < TOL)) {
      which_0 <- gf < TOL
      for (j in seq_len(nrow(xmat))) {
        if (!all(xmat[j, which_0] == 0)) {
          expect_true(is_impossible(x = xmat[j, ], g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], dr = pdf$dr[[i]]))
        } else {
          expect_false(is_impossible(x = xmat[j, ], g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], dr = pdf$dr[[i]]))
        }
      }
    } else {
      expect_false(is_impossible(x = xmat[j, ], g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], dr = pdf$dr[[i]]))
    }
  }
})


test_that("is_valid_2() works", {
  expect_true(is_valid_2(dr = 1/6, pp = 1/3, drbound = 1/6))
  expect_false(is_valid_2(dr = 1/6, pp = 1/3 - 1e-6, drbound = 1/6))
  expect_false(is_valid_2(dr = 1/6, pp = 1/3 + 1e-6, drbound = 1/6))
  expect_true(is_valid_2(dr = 0, pp = 0, drbound = 1/6))
  expect_true(is_valid_2(dr = 0, pp = 1, drbound = 1/6))
})

test_that("three and two equal", {
  tau <- 1/2
  beta <- 1/6
  gamma <- 1/3
  ell <- 2
  p3 <- pvec_tet_3(tau = tau, beta = beta, gamma = gamma, ell = ell)
  par <- three_to_two(tau = tau, beta = beta, gamma = gamma)
  alpha <- par[[1]]
  xi <- par[[2]]
  p2 <- pvec_tet_2(alpha = alpha, xi = xi, ell = ell)
  expect_equal(p2, p3)
})




