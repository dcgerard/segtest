test_that("null simulations produce uniform p-values", {
  skip("not for unit testing")
  p1_ploidy <- 4
  p1 <- 1
  p2_ploidy <- 4
  p2 <- 1
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
    pi = 0.01)
  niter <- 1000
  nsamp <- 10000


  pval_vec <- rep(NA_real_, times = niter)
  df_vec <- rep(NA_real_, times = niter)
  stat_vec <- rep(NA_real_, times = niter)
  df1_vec <- rep(NA_real_, times = niter)
  df0_vec <- rep(NA_real_, times = niter)
  pi_vec <- rep(NA_real_, times = niter)

  for (i in seq_len(niter)) {
    if (i %% 50 == 0) {
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



