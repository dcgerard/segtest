test_that("hard cases work", {
  x <- c(6L, 87L, 429L, 392L, 86L)
  g1 <- 3
  g2 <- 2
  par = c(1e-5, 1e-5, 1 - 1e-5)
  expect_no_error(
    obj_dr_pp(par = par, x = x, g1 = g1, g2 = g2)
  )
})

test_that("hard case works for tetraploids", {
  # load("./mout_Cc.rda")
  # mout_Cc$snpdf$snp[3055:3056]
  # sub <- filter_snp(mout_Cc, snp %in% c("chr4c_2876569", "chr4c_2876589"))
  # mhard_g <- multidog_to_g(mout = sub, ploidy = 4, type = "off_g")
  # mhard_gl <- multidog_to_g(mout = sub, ploidy = 4, type = "off_gl")
  # save(mhard_g, mhard_gl, file = "./mhard.RData")

  load("./mhard.RData")

  expect_no_error(
    {
      #seg_lrt(x = mhard_gl$g[1, ,], p1_ploidy = 4, p2_ploidy = 4, p1 = mhard_gl$p1[[1]], p2 = mhard_gl$p1[[1]])
      #seg_multi(g = mhard_gl$g, p1_ploidy = 4, p2_ploidy = 4, p1 = mhard_gl$p1, p2 = mhard_gl$p2)
      em_li(mhard_gl$g[1, ,])
    }
  )
})
