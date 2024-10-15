test_that("multi_lrt gives same stuff as lrt funs", {
  i <- 3
  ## Assuming genotypes are known (typically bad idea)
  glist <- multidog_to_g(mout = ufit, type = "all_g", p1 = "indigocrisp", p2 = "sweetcrisp")
  p1_1 <- glist$p1
  p2_1 <- glist$p2
  g_1 <- glist$g
  mout <- multi_lrt(g = g_1, p1 = p1_1, p2 = p2_1)
  uout <- unlist(lrt_men_g4(x = g_1[i, ], g1 = p1_1[[i]], g2 = p2_1[[i]]))
  expect_equal(mout$p_value[[i]], uout[["p_value"]])

  ## Using genotype likelihoods (typically good idea)
  glist <- multidog_to_g(mout = ufit, type = "all_gl", p1 = "indigocrisp", p2 = "sweetcrisp")
  p1_2 <- glist$p1
  p2_2 <- glist$p2
  g_2 <- glist$g
  mout <- multi_lrt(g = g_2, p1 = p1_2, p2 = p2_2)
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


test_that("Get uniform distribution under null with multi_lrt()", {
  skip("too long")

  ## Change these
  g1 <- 2
  g2 <- 2
  dr <- TRUE
  pp <- TRUE
  rd <- Inf
  n <- 100
  alpha <- 1/12
  xi1 <- 0
  xi2 <- 1/4
  nloc <- 1000

  snpname <- paste0("Ind", seq_len(nloc))

  if (is.infinite(rd)) {
    g <- replicate(n = nloc, expr = {
      simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    })
    g <- t(g)
    rownames(g) <- snpname
    p1 <- rep(g1, nloc)
    names(p1) <- snpname
    p2 <- rep(g1, nloc)
    names(p2) <- snpname
  } else {
    g <- replicate(n = nloc, expr = {
      suppressWarnings(
        simf1gl(n = n, g1 = g1, g2 = g2, rd = rd, alpha = alpha, xi1 = xi1, xi2 = xi2)
      )
    })


  }

})
