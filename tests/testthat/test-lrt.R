test_that("LRT works on edge cases", {
  x <- c(10, 0, 0, 0, 0)

  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = TRUE, dr = FALSE)$p_value != 0)
  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = FALSE, dr = FALSE)$p_value != 0)
})


test_that("qq plot is unif in some cases", {
  skip("too slow")


  g1 <- 0
  g2 <- 2
  alpha <- 0
  xi1 <- 0
  xi2 <- 1
  pp <- TRUE
  dr <- TRUE
  n <- 1000
  iter <- 100
  pvec <- rep(NA_real_, iter)
  stat <- rep(NA_real_, iter)
  df <- rep(NA_real_, iter)
  aest <- rep(NA_real_, iter)
  xi1est <- rep(NA_real_, iter)
  xi2est <- rep(NA_real_, iter)
  gf <- matrix(ncol = 5, nrow = iter)
  for (i in seq_len(iter)) {
    if (i %% 10 == 0) {
      cat(i, " of ", iter, "\n")
    }
    x <- simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    lout <- lrt_men_g4(x = x, g1 = g1, g2 = g2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
    pvec[[i]] <- lout$p_value
    stat[[i]] <- lout$statistic
    df[[i]] <- lout$df
    #gf[i, ] <- offspring_gf_2(alpha = lout$alpha, xi1 = lout$xi1, xi2 = lout$xi2, p1 = g1, p2 = g2)
    aest[[i]] <- lout$alpha
    xi1est[[i]] <- lout$xi1
    xi2est[[i]] <- lout$xi2
  }

  table(df)

  ## hwep::qqpvalue(pvals = pvec)
  qqplot(x = ppoints(iter), y = pvec, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  ## competitors
  pvec2 <- stats::pchisq(q = stat, df = 2, lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec2, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  pvec3 <- stats::pchisq(q = stat, df = ifelse(aest < alpha - 1e-7, 4, 3), lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec3, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)
})


test_that("GL qq plot is unif in some cases", {
  skip("too slow")


  g1 <- 0
  g2 <- 4
  alpha <- 0
  xi1 <- 0
  xi2 <- 0
  pp <- TRUE
  dr <- TRUE
  rd <- 10
  n <- 200
  iter <- 100
  pvec <- rep(NA_real_, iter)
  stat <- rep(NA_real_, iter)
  df <- rep(NA_real_, iter)
  aest <- rep(NA_real_, iter)
  xi1est <- rep(NA_real_, iter)
  xi2est <- rep(NA_real_, iter)
  gf <- matrix(ncol = 5, nrow = iter)
  for (i in seq_len(iter)) {
    cat(i, " of ", iter, "\n")
    gl <- simf1gl(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, rd = rd)
    lout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
    pvec[[i]] <- lout$p_value
    stat[[i]] <- lout$statistic
    df[[i]] <- lout$df
    #gf[i, ] <- offspring_gf_2(alpha = lout$alpha, xi1 = lout$xi1, xi2 = lout$xi2, p1 = g1, p2 = g2)
    aest[[i]] <- lout$alpha
    xi1est[[i]] <- lout$xi1
    xi2est[[i]] <- lout$xi2
  }

  table(df)

  ## hwep::qqpvalue(pvals = pvec)
  qqplot(x = ppoints(iter), y = pvec, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  ## competitors
  pvec2 <- stats::pchisq(q = stat, df = 3, lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec2, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  pvec3 <- stats::pchisq(q = stat, df = ifelse(aest > 1/6 - 1e-7, 2, 1), lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec3, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)
})


test_that("corner cases", {
  if (!requireNamespace("polymapR", quietly = TRUE)) {
    skip("polymapR not installed")
  }
  set.seed(1)

  xmat <- structure(c(4L, 8L, 8L, 9L, 7L, 18L, 8L, 7L, 10L, 1L, 0L, 0L,
0L, 1L, 1L, 2L, 3L, 4L, 4L, 5L, 3L, 5L, 2L, 6L, 4L, 5L, 3L, 2L,
2L, 3L, 6L, 5L, 4L, 7L, 10L, 7L, 8L, 8L, 0L, 0L, 7L, 2L, 0L,
0L, 0L, 0L, 0L, 5L, 0L, 0L, 1L, 1L, 0L, 3L, 3L, 2L, 8L, 16L,
10L, 9L, 8L, 11L, 1L, 8L, 10L, 7L, 8L, 19L, 18L, 20L, 18L, 19L,
18L, 17L, 5L, 7L, 9L, 9L, 8L, 13L, 6L, 8L, 10L, 12L, 9L, 12L,
11L, 7L, 6L, 9L, 8L, 8L, 5L, 6L, 7L, 4L, 16L, 9L, 6L, 3L, 2L,
0L, 1L, 0L, 6L, 0L, 1L, 1L, 2L, 1L, 56L, 44L, 53L, 36L, 0L, 2L,
3L, 3L, 2L, 1L, 4L, 3L, 3L, 11L, 1L, 2L, 0L, 1L, 0L, 0L, 0L,
11L, 7L, 4L, 5L, 5L, 3L, 6L, 5L, 3L, 3L, 6L, 4L, 4L, 5L, 7L,
5L, 2L, 0L, 6L, 5L, 4L, 16L, 4L, 2L, 9L, 3L, 15L, 18L, 18L, 20L,
6L, 16L, 15L, 15L, 15L, 16L, 94L, 87L, 97L, 121L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L,
2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 2L, 2L, 2L, 2L, 1L, 3L, 2L,
2L, 1L, 1L, 0L, 0L, 2L, 1L, 12L, 2L, 2L, 1L, 0L, 3L, 3L, 3L,
3L, 2L, 2L, 44L, 63L, 44L, 32L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 2L, 2L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 3L, 3L, 4L,
3L), dim = c(57L, 5L), dimnames = list(NULL, c("0", "1", "2",
"3", "4")))

  ell1vec <- c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
  ell2vec <- c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)

  for(i in seq_len(nrow(xmat))) {
    x <- xmat[i, ]
    g1 <- ell1vec[[i]]
    g2 <- ell2vec[[i]]
    mout <- polymapr_test(x = x, g1 = g1, g2 = g2, type = "segtest")$p_value
    pout <- polymapr_test(x = x, g1 = g1, g2 = g2, type = "polymapR")$p_value
    expect_equal(mout - pout, 0, tolerance = 10^-3)
  }

  ## on last i
  gl <- simgl(nvec = x)
  pl <- t(apply(X = gl, MARGIN = 1, FUN = \(x) exp(x - log_sum_exp(x))))
  mout <- polymapr_test(x = pl, g1 = g1, g2 = g2, type = "segtest")$p_value
  pout <- polymapr_test(x = pl, g1 = g1, g2 = g2, type = "polymapR")$p_value
  expect_equal(mout - pout, 0, tolerance = 10^-3)
})

