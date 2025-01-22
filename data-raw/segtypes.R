## code to prepare `segtypes` dataset goes here
library(polymapR)
library(tidyverse)
col <- polymapR::calcSegtypeInfo(ploidy = 4, ploidy2 = 4)
tibble(col = col) |>
  unnest_wider(col = col) |>
  mutate(mod = names(col)) ->
  segtypes

usethis::use_data(segtypes, overwrite = TRUE, internal = TRUE)


## code to prepare `seg` dataset goes here
slist <- list()
for (i in 1:10) {
  ploidy <- i * 2
  df_allo <- pp_seg_pats(ploidy = ploidy)
  df_auto <- data.frame(
    ploidy = ploidy,
    g = 0:ploidy,
    m = NA,
    p = I(vector(mode = "list", length = ploidy + 1)))
  for (ell in 0:ploidy) {
    df_auto$p[[ell + 1]] <- list()
    df_auto$p[[ell + 1]] <- stats::dhyper(x = 0:(ploidy / 2), m = ell, n = ploidy - ell, k = ploidy / 2)
  }
  df_allo$mode <- "disomic"
  df_allo$mode[df_allo$g == 0 | df_allo$g == 1 | df_allo$g == ploidy - 1 | df_allo$g == ploidy] <- "both"
  df_auto$mode <- "polysomic"

  df <- rbind(df_allo, df_auto)
  df <- df[!(df$mode == "polysomic" & (df$g == 0 | df$g == 1 | df$g == ploidy - 1 | df$g == ploidy)), ]
  slist[[i]] <- df
}

seg <- do.call(what = rbind, args = slist)
usethis::use_data(seg)




