test_that("null simulations produce uniform p-values", {
  skip("not for unit testing")
  p1_ploidy <- 4
  p1 <- 1
  p2_ploidy <- 4
  p2 <- 3
  q <- gf_freq(
    p1_g = p1,
    p1_ploidy = p1_ploidy,
    p1_gamma = 1,
    p1_beta = 0.02,
    p1_alpha = NULL,
    p1_type = "mix",
    p2_g = p2,
    p2_ploidy = p2_ploidy,
    p2_gamma= 1,
    p2_beta = 0.02,
    p2_alpha = NULL,
    p2_type = "mix",
    pi = 0.015)
  niter <- 1000
  nsamp <- 10000


  pval_vec <- rep(NA_real_, times = niter)
  df_vec <- rep(NA_real_, times = niter)
  stat_vec <- rep(NA_real_, times = niter)
  df1_vec <- rep(NA_real_, times = niter)
  df0_vec <- rep(NA_real_, times = niter)
  pi_vec <- rep(NA_real_, times = niter)

  for (i in seq_len(niter)) {
    if (i %% 10 == 0) {
      cat(i, "out of", niter, "\n")
    }
    nvec <- c(stats::rmultinom(n = 1, size = nsamp, prob = q))
    sout <- seg_lrt(
      x = nvec,
      p1_ploidy = p1_ploidy,
      p2_ploidy = p2_ploidy,
      p1 = p1,
      p2 = p2,
      model = "seg",
      outlier = TRUE)
    pval_vec[[i]] <- sout$p_value
    df_vec[[i]] <- sout$df
    stat_vec[[i]] <- sout$stat
    df1_vec[[i]] <- sout$alt$df1
    df0_vec[[i]] <- sout$null$df0
    pi_vec[[i]] <- sout$null$gam[[3]]$pi
  }

  plot(
    x = ppoints(n = sum(!is.na(pval_vec))),
    y = sort(pval_vec[!is.na(pval_vec)]),
    xlab = "Theoretical",
    ylab = "Sample")
  abline(a = 0, b = 1, lty = 2, col = 2)
})

test_that("corner cases work", {
  sout <- seg_lrt(x = c(1, 0, 0, 0, 0), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 0, outlier = FALSE)
  expect_equal(sout$p_value, 1, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 0, 1, 0, 0), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 4, outlier = FALSE)
  expect_equal(sout$p_value, 1, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 0, 0, 0, 1), p1_ploidy = 4, p2_ploidy = 4, p1 = 4, p2 = 4, outlier = FALSE)
  expect_equal(sout$p_value, 1, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 1, 0, 0, 0), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 0, outlier = FALSE)
  expect_equal(sout$p_value, 0, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 0, 1, 0, 0), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 0, outlier = FALSE)
  expect_equal(sout$p_value, 0, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 0, 0, 1, 0), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 0, outlier = FALSE)
  expect_equal(sout$p_value, 0, tolerance = 1e-3)

  sout <- seg_lrt(x = c(0, 0, 0, 0, 1), p1_ploidy = 4, p2_ploidy = 4, p1 = 0, p2 = 0, outlier = FALSE)
  expect_equal(sout$p_value, 0, tolerance = 1e-3)
})

test_that("skip", {
  skip("testing out frequencies")
  library(numDeriv)
  fn <- function(par) {
    gf_freq(
      p1_g = 3,
      p1_ploidy = 4,
      p1_gamma = 1,
      p1_beta = par[[1]],
      p1_alpha = NULL,
      p1_type = "mix",
      p2_g = 1,
      p2_ploidy = 4,
      p2_gamma= 1,
      p2_beta = par[[2]],
      p2_alpha = NULL,
      p2_type = "mix",
      pi = par[[3]])
  }
  q <- fn(par = c(0, 0, 0.015))
  par <- c(beta_bounds(4), beta_bounds(4), 0.03)
  svd(numDeriv::jacobian(func = fn, x = par/2))$d

  obj <- function(par) {
    sum((fn(par) - q)^2)
  }
  oout <- optim(par = par, fn = obj, lower = c(0, 0, 0), upper = c(0.1, 0.1, 0.01), method = "L-BFGS-B")

  oout$par
  fn(oout$par) - q
  q
})

test_that("seg at simplex and auto_dr are same", {
  q1 <- gamfreq(g = 1, ploidy = 4, gamma = 1, alpha = NULL, beta = 0.01, type = "mix", add_dr = TRUE)
  q2 <- gamfreq(g = 1, ploidy = 4, gamma = NULL, alpha = 0.04, beta = NULL, type = "polysomic", add_dr = TRUE)
  expect_equal(q1, q2)

  q1 <- gamfreq(g = 3, ploidy = 4, gamma = 1, alpha = NULL, beta = 0.01, type = "mix", add_dr = TRUE)
  q2 <- gamfreq(g = 3, ploidy = 4, gamma = NULL, alpha = 0.04, beta = NULL, type = "polysomic", add_dr = TRUE)
  expect_equal(q1, q2)

  q1 <- gf_freq(p1_g = 1, p1_ploidy = 4, p1_gamma = 1, p1_alpha = NULL, p1_beta = 0.01, p1_type = "mix", p1_add_dr = TRUE, p2_g = 3, p2_ploidy = 4, p2_gamma = 1, p2_alpha = NULL, p2_beta = 0.01, p2_type = "mix", p2_add_dr = TRUE, pi = 0.01)
  q2 <- gf_freq(p1_g = 1, p1_ploidy = 4, p1_gamma = NULL, p1_alpha = 0.04, p1_beta = NULL, p1_type = "polysomic", p1_add_dr = TRUE, p2_g = 3, p2_ploidy = 4, p2_gamma = NULL, p2_alpha = 0.04, p2_beta = NULL, p2_type = "polysomic", p2_add_dr = TRUE, pi = 0.01)
  expect_equal(q1, q2)

  q1 <- gf_freq(p1_g = 3, p1_ploidy = 4, p1_gamma = 1, p1_alpha = NULL, p1_beta = 0.01, p1_type = "mix", p1_add_dr = TRUE, p2_g = 1, p2_ploidy = 4, p2_gamma = 1, p2_alpha = NULL, p2_beta = 0.01, p2_type = "mix", p2_add_dr = TRUE, pi = 0.01)
  q2 <- gf_freq(p1_g = 3, p1_ploidy = 4, p1_gamma = NULL, p1_alpha = 0.04, p1_beta = NULL, p1_type = "polysomic", p1_add_dr = TRUE, p2_g = 1, p2_ploidy = 4, p2_gamma = NULL, p2_alpha = 0.04, p2_beta = NULL, p2_type = "polysomic", p2_add_dr = TRUE, pi = 0.01)
  expect_equal(q1, q2)
})

test_that("auto_dr and seg are same", {
  set.seed(1)
  nvec <- c(32L, 2431L, 4979L, 2520L, 38L)
  sout1 <- seg_lrt(
    x = nvec,
    p1_ploidy = 4,
    p2_ploidy = 4,
    p1 = 3,
    p2 = 1,
    model = "auto_dr",
    outlier = TRUE)
  sout2 <- seg_lrt(
    x = nvec,
    p1_ploidy = 4,
    p2_ploidy = 4,
    p1 = 3,
    p2 = 1,
    model = "seg",
    outlier = TRUE)
  expect_equal(sout1$null$l0, sout2$null$l0)
})
