## code to prepare `segtypes` dataset goes here
library(polymapR)
library(tidyverse)
col <- polymapR::calcSegtypeInfo(ploidy = 4, ploidy2 = 4)
tibble(col = col) |>
  unnest_wider(col = col) |>
  mutate(mod = names(col)) ->
  segtypes

usethis::use_data(segtypes, overwrite = TRUE, internal = TRUE)
