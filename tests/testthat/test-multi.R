test_that("multi_lrt gives same stuff as lrt funs", {
  i <- 3
  ## Assuming genotypes are known (typically bad idea)
  glist <- multidog_to_g(mout = ufit, ploidy = 4, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
  p1_1 <- glist$p1
  p2_1 <- glist$p2
  g_1 <- glist$g
  mout <- multi_lrt(g = g_1, p1 = p1_1, p2 = p2_1)
  uout <- unlist(lrt_men_g4(x = g_1[i, ], g1 = p1_1[[i]], g2 = p2_1[[i]]))
  expect_equal(mout$p_value[[i]], uout[["p_value"]])

  ## Using genotype likelihoods (typically good idea)
  glist <- multidog_to_g(mout = ufit, ploidy = 4, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
  p1_2 <- glist$p1
  p2_2 <- glist$p2
  g_2 <- glist$g
  mout <- multi_lrt(g = g_2, p1 = p1_2, p2 = p2_2, nullprop = TRUE)
  uout <- unlist(lrt_men_gl4(gl = g_2[i,,], g1 = p1_2[i,], g2 = p2_2[i,]))
  expect_equal(mout$p_value[[i]], uout[["p_value"]])

  ## Offspring genotype likelihoods and parent genotypes known
  mout <- multi_lrt(g = g_2, p1 = p1_1, p2 = p2_1)
  uout <- unlist(lrt_men_gl4(gl = g_2[i,,], g1 = p1_1[[i]], g2 = p2_2[[i]]))
  expect_equal(mout$p_value[[i]], uout[["p_value"]])

  ## Offspring genotype likelihoods and no information on parent genotypes
  mout <- multi_lrt(g = g_2, p1 = NULL, p2 = NULL)
  uout <- unlist(lrt_men_gl4(gl = g_2[i,,], g1 = NULL, g2 = NULL))
  expect_equal(mout$p_value[[i]], uout[["p_value"]])
})

test_that("seg_multi() gives same stuff as seg_lrt()", {
  set.seed(1)
  i <- 5
  glist <- multidog_to_g(mout = ufit, ploidy = 4, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
  p1_1 <- glist$p1
  p2_1 <- glist$p2
  g_1 <- glist$g
  mout <- seg_multi(g = g_1, p1_ploidy = 4, p2_ploidy = 4, p1 = p1_1, p2 = p2_1)
  sout <- seg_lrt(x = g_1[i, ], p1_ploidy = 4, p2_ploidy = 4, p1 = p1_1[[i]], p2 = p2_1[[i]])
  expect_equal(
    mout$p_value[[i]],
    sout$p_value,
    tolerance = 1e-5
  )
})


test_that("Get uniform distribution under null with multi_lrt()", {
  skip("too long")

  ## Change these
  g1 <- 2
  g2 <- 1
  dr <- TRUE
  pp <- FALSE
  rd <- Inf
  n <- 1000
  alpha <- 1/12
  xi1 <- 1/3
  xi2 <- 1/3
  nloc <- 1000

  is_valid_2(dr = alpha, pp = xi1)

  snpname <- paste0("Loc", seq_len(nloc))

  if (is.infinite(rd)) {
    g <- replicate(n = nloc, expr = {
      simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    })
    g <- t(g)
    rownames(g) <- snpname
    p1 <- rep(g1, nloc)
    names(p1) <- snpname
    p2 <- rep(g2, nloc)
    names(p2) <- snpname
  } else {
    g <- replicate(n = nloc, expr = {
      suppressWarnings(
        simf1gl(n = n, g1 = g1, g2 = g2, rd = rd, alpha = alpha, xi1 = xi1, xi2 = xi2)
      )
    })
    g <- aperm(a = g, perm = c(3, 1, 2))
    indname <- paste0("Ind", seq_len(n))
    dimnames(g) <- list(snpname, indname, 0:4)
    p1 <- rep(g1, nloc)
    names(p1) <- snpname
    p2 <- rep(g1, nloc)
    names(p2) <- snpname
  }

  future::plan(future::multisession(workers = 2))
  mout <- multi_lrt(g = g, p1 = p1, p2 = p2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
  future::plan(future::sequential())

  stats::qqplot(
    x = stats::ppoints(nloc),
    y = mout$p_value,
    main = "QQ Plot of P-values",
    xlab = "Theoretical Quantiles",
    ylab = "P-values")
  graphics::text(x = 0.6, y = 0.1, label = "Valid if at or above y=x line", col = "red")
  graphics::abline(a = 0, b = 1)

})


test_that("multidog_to_g works with S1 data", {
  load("./s1dat.RData")
  ## uout_s1_1 fit with model s1 but no parent given
  ## uout_s1_2 fit with model s1 but a parent was provided
  ## uout_s1_3 fit with model norm, parent id is Xushu18

  expect_no_error(
    {
      multidog_to_g(mout = uout_s1_1, ploidy = 6, type = "off_g")
      multidog_to_g(mout = uout_s1_1, ploidy = 6, type = "off_gl")
      multidog_to_g(mout = uout_s1_2, ploidy = 6, type = "off_g")
      multidog_to_g(mout = uout_s1_2, ploidy = 6, type = "off_gl")
      multidog_to_g(mout = uout_s1_3, ploidy = 6, type = "all_g", p1 = "Xushu18")
      multidog_to_g(mout = uout_s1_3, ploidy = 6, type = "all_gl", p1 = "Xushu18")
    })
})


test_that("multidog_to_g works with ploidy = 4 for any model", {
  load("./all4_dat.RData")

  expect_no_error(
    {
      multidog_to_g(ufit_s1, ploidy = 4, type = "off_g")
      multidog_to_g(ufit_s1, ploidy = 4, type = "off_gl")
      multidog_to_g(ufit_s1pp, ploidy = 4, type = "off_g")
      multidog_to_g(ufit_s1pp, ploidy = 4, type = "off_gl")
      multidog_to_g(ufit_f1, ploidy = 4, type = "off_g")
      multidog_to_g(ufit_f1, ploidy = 4, type = "off_gl")
      multidog_to_g(ufit_f1pp, ploidy = 4, type = "off_g")
      multidog_to_g(ufit_f1pp, ploidy = 4, type = "off_gl")
      multidog_to_g(ufit_norm, ploidy = 4, type = "all_g", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_norm, ploidy = 4, type = "all_g", p1 = "sweetcrisp")
      multidog_to_g(ufit_norm, ploidy = 4, type = "all_gl", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_norm, ploidy = 4, type = "all_gl", p1 = "sweetcrisp")
      multidog_to_g(ufit_hw, ploidy = 4, type = "all_g", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_hw, ploidy = 4, type = "all_g", p1 = "sweetcrisp")
      multidog_to_g(ufit_hw, ploidy = 4, type = "all_gl", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_hw, ploidy = 4, type = "all_gl", p1 = "sweetcrisp")
      multidog_to_g(ufit_bb, ploidy = 4, type = "all_g", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_bb, ploidy = 4, type = "all_g", p1 = "sweetcrisp")
      multidog_to_g(ufit_bb, ploidy = 4, type = "all_gl", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_bb, ploidy = 4, type = "all_gl", p1 = "sweetcrisp")
      multidog_to_g(ufit_flex, ploidy = 4, type = "all_g", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_flex, ploidy = 4, type = "all_g", p1 = "sweetcrisp")
      multidog_to_g(ufit_flex, ploidy = 4, type = "all_gl", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_flex, ploidy = 4, type = "all_gl", p1 = "sweetcrisp")
      multidog_to_g(ufit_uniform, ploidy = 4, type = "all_g", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_uniform, ploidy = 4, type = "all_g", p1 = "sweetcrisp")
      multidog_to_g(ufit_uniform, ploidy = 4, type = "all_gl", p1 = "sweetcrisp", p2 = "indigocrisp")
      multidog_to_g(ufit_uniform, ploidy = 4, type = "all_gl", p1 = "sweetcrisp")
    })

})


