test_that("which_possible() works", {
  df <- expand.grid(g1 = 0:4, g2 = 0:4)
  for (i in seq_len(nrow(df))) {
    expect_false(is_impossible(which_possible(g1 = df$g1[[i]], g2 = df$g2[[i]], dr = TRUE) * 1, g1 = df$g1[[i]], g2 = df$g2[[i]], dr = TRUE))
    expect_false(is_impossible(which_possible(g1 = df$g1[[i]], g2 = df$g2[[i]], dr = FALSE) * 1, g1 = df$g1[[i]], g2 = df$g2[[i]], dr = FALSE))
  }
})
