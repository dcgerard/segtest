test_that("drbounds works", {
  expect_equal(drbounds(ploidy = 4, model = "ces"), 1/6)
  expect_equal(drbounds(ploidy = 6, model = "ces"), 3/10)
  expect_equal(drbounds(ploidy = 8, model = "ces"), c(54, 3) / 140)
  expect_equal(drbounds(ploidy = 10, model = "ces"), c(110, 15) / 252)
  expect_equal(drbounds(ploidy = 12, model = "ces"), c(855, 195, 5) / 1848)

  expect_equal(drbounds(ploidy = 4, model = "prcs"), 1/7)
  expect_equal(drbounds(ploidy = 6, model = "prcs"), 3/11)
  expect_equal(drbounds(ploidy = 8, model = "prcs"), c(24, 1) / 65)
  expect_equal(drbounds(ploidy = 10, model = "prcs"), c(140, 15) / 323)
  expect_equal(drbounds(ploidy = 12, model = "prcs"), c(1440, 270, 5) / 3059)

  expect_equal(beta_bounds(ploidy = 4, model = "ces"), expected = 1/24)
  expect_equal(beta_bounds(ploidy = 6, model = "ces"), expected = 1/20)
  expect_equal(beta_bounds(ploidy = 8, model = "ces"), expected = 3/56)
  expect_equal(beta_bounds(ploidy = 10, model = "ces"), expected = 1/18)
  expect_equal(beta_bounds(ploidy = 12, model = "ces"), expected = 5/88)
})
