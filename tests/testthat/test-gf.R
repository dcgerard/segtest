test_that("old and new gf are the same", {

  for (g1 in 0:4) {
    for (g2 in 0:4) {
      gf1 <- gf_freq(
        p1_g = g1,
        p1_ploidy = ploidy,
        p1_gamma = NULL,
        p1_alpha = 0,
        p1_beta = NULL,
        p1_type = "polysomic",
        p1_add_dr = FALSE,
        p2_g = g2,
        p2_ploidy = ploidy,
        p2_gamma = NULL,
        p2_alpha = 0,
        p2_beta = NULL,
        p2_type = "polysomic",
        p2_add_dr = FALSE,
        pi = 0)
      gf2 <- offspring_gf_2(
        alpha = 0,
        xi1 = 1/3,
        xi2 = 1/3,
        p1 = g1,
        p2 = g2)
      expect_equal(gf1, gf2, tolerance = 1e-6)
    }
  }
})
