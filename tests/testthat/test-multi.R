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
