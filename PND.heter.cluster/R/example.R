
# library(tidyverse)
# library(SuperLearner)
# library(glue)
#
# data(data_in)
# data_in <- data_in
# Xnames <- c(grep("X_dat", colnames(data_in), value = TRUE))
#
# estimates_ate_K <- atekCl(
# data_in = data_in, ttname = "tt", Kname = "K", Yname = "Y",
# Xnames = Xnames,
# Yfamily = "gaussian",
# sensitivity = NULL # or "small_to_medium"
# )
# estimates_ate_K$ate_K
