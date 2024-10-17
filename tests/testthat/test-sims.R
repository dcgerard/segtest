test_that("gl sims work", {
  set.seed(1)
  sout <- simf1gl(n = 10000, g1 = 2, g2 = 2, rd = 1000, alpha = 1/6, xi1 = 1/3, xi2 = 1/3)
  gf_true <- offspring_gf_2(alpha = 1/6, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
  x <- table(factor(apply(sout, 1, which.max) - 1, levels = 0:4))
  gf_emp <- x / sum(x)
  expect_true(sum(abs(gf_true - gf_emp)) < 0.1)
})

test_that("offpsring_geno works", {
  gf_true <- offspring_gf_2(alpha = 1/6, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
  x <- offspring_geno(gf = gf_true, n = 10000)
  gf_emp <- x / sum(x)
  expect_true(sum(abs(gf_true - gf_emp)) < 0.1)
})
