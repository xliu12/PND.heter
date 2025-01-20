# PND.heter.cluster
R package for assessing heterogeneous cluster-specific treatment effects in partially nested designs (also known as "partially clustered designs"). 

For details, please see the accompanying methods paper:

Liu, X. (2024). Assessing heterogeneous causal effects across clusters in partially nested designs. Psychological Methods. Advance online publication. https://doi.org/10.1037/met0000723

## Installation

To install `PND.heter.cluster` directly from Github:
```
remotes::install_github(repo = "xliu12/PND.heter", subdir = "PND.heter.cluster")
```

## Example

```
library(tidyverse)
library(SuperLearner)
library(glue)

# data
data(data_in)
data_in <- data_in

# baseline covariates
Xnames <- c(grep("X_dat", colnames(data_in), value = TRUE))

estimates_ate_K <- atekCl(
  data_in = data_in,
  ttname = "tt",  # treatment variable
  Kname = "K",    # cluster assignment variable, coded as 0 for individuals in the (non-clustered) control arm
  Yname = "Y",    # outcome variable
  Xnames = Xnames
)
estimates_ate_K$ate_K

```
