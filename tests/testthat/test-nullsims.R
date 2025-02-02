test_that("null simulations produce uniform p-values", {
  skip("not for unit testing")
  p1_ploidy <- 6
  p1 <- 3
  p2_ploidy <- 6
  p2 <- 3
  q <- gf_freq(
    p1_g = p1,
    p1_ploidy = p1_ploidy,
    p1_gamma = c(0.5, 0.5),
    p1_beta = NULL,
    p1_type = "mix",
    p2_g = p2,
    p2_ploidy = p2_ploidy,
    p2_gamma= c(0.3, 0.7),
    p2_type = "mix",
    pi = 0.03)
  niter <- 1000
  nsamp <- 200


  pval_vec <- rep(NA_real_, times = niter)
  df_vec <- rep(NA_real_, times = niter)
  stat_vec <- rep(NA_real_, times = niter)

  for (i in seq_len(niter)) {
    if (i %% 50 == 0) {
      cat(i, "out of", niter, "\n")
    }
    nvec <- c(stats::rmultinom(n = 1, size = nsamp, prob = q))
    sout <- seg_lrt(x = nvec, p1_ploidy = p1_ploidy, p2_ploidy = p2_ploidy, p1 = p1, p2 = p2)
    pval_vec[[i]] <- sout$p_value
    df_vec[[i]] <- sout$df
    stat_vec[[i]] <- sout$stat
  }

  plot(
    x = ppoints(n = sum(!is.na(pval_vec))),
    y = sort(pval_vec[!is.na(pval_vec)]),
    xlab = "Theoretical",
    ylab = "Sample")
  abline(a = 0, b = 1, lty = 2, col = 2)


})
