test_that("Corner cases work", {
  x <- c(`0` = 95L, `1` = 397L, `2` = 417L, `3` = 87L, `4` = 4L)
  g1 <- 2
  g2 <- 1
  alpha <- 1/12
  xi1 <- 1/3
  xi2 <- 1/3
  pp <- TRUE
  dr <- TRUE
  expect_no_error({
  lout <- lrt_men_g4(
    x = x,
    g1 = g1,
    g2 = g2,
    pp = pp,
    dr = dr,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)
  })

  obj_dr_pp(par = c(0, 0, 0.999), x = x, g1 = g1, g2 = g2)

  like_gknown_3(
      x = x,
      tau = 1e-7,
      beta = 1e-7,
      gamma1 = 0.999,
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
})

test_that("corner cases work for higher ploidies", {
  expect_no_error(
    seg_lrt(x = c(1, 0, 0, 0, 0), p1_ploidy = 4, p2_ploidy = 4)
  )

  expect_no_error(
    seg_lrt(x = c(1, 0, 0, 1, 0), p1_ploidy = 4, p2_ploidy = 4)
  )

  expect_no_error(
    seg_lrt(x = c(0, 1, 0, 1, 0, 0, 0), p1_ploidy = 6, p2_ploidy = 6)
  )
})


test_that("new corner case works", {
  load("./fpop.RData")
  ploidy <- 6
  p1geno <- fpop$par$p1geno
  p2geno <- fpop$par$p2geno
  gl <- fpop$genologlike
  sout <- seg_lrt(x = gl, p1_ploidy = ploidy, p2_ploidy = ploidy, p1 = p1geno, p2 = p2geno)
  expect_true(!is.na(sout$p_value))
})



