test_that("lrt_men_g4 works", {
  expect_no_error({
    pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
    for (i in seq_len(nrow(pdf))) {
      set.seed(i)
      x <- round(100 * offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = pdf$p1[[i]], p2 = pdf$p2[[i]]))
      t1 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE)
      t2 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE)
      t3 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE)
      t4 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE)
      t5 <- chisq_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]])
    }
  })
})

test_that("lrt_men_gl4 works", {
  expect_no_error({
    pdf <- data.frame(p1 = sample(0:4, size = 1), p2 = sample(0:4, size = 1))
    for (i in seq_len(nrow(pdf))) {
      set.seed(i)
      gl <- simf1gl(n = 10, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]])
      t1 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE)
      t2 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE)
      t3 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE)
      t4 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE)
      t5 <- lrt_men_gl4(gl, pp = TRUE, dr = TRUE)
      t6 <- lrt_men_gl4(gl, pp = TRUE, dr = FALSE)
      t7 <- lrt_men_gl4(gl, pp = FALSE, dr = TRUE)
      t8 <- lrt_men_gl4(gl, pp = FALSE, dr = FALSE)
      t9 <- chisq_gl4(gl = gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]])
    }
  })
})

test_that("otest_g gives same as lrt", {
  x <- c(5, 26, 40, 24, 10)
  expect_equal(otest_g(x = x, g1 = 2, g2 = 2)$p_lrt, lrt_men_g4(x = x, g1 = 2, g2 = 2)$p_value)
})

