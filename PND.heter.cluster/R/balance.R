
#' Checking covariate balance based on estimated cluster assignment probabilities (principal score) and treatment assignment probabilities (propensity score).
#' @param data_in A \code{data.frame} containing all necessary variables.
#' @param atekCl_results [\code{list}]\cr
#'  A list returned from the R function \code{atekCl()}.
#' @param Kname [\code{character}]\cr
#'  A character string of the column name of the cluster assignment variable. This variable should be coded as 0 for individuals in the control arm, the arm without the cluster assignment.
#' @param ttname [\code{character}]\cr
#'  A character string of the column name of the treatment variable. The treatment variable should be dummy-coded, with 1 for the (clustered) treatment arm and 0 for the (non-clustered) control arm.
#' @param covariate_names  [\code{character}]\cr
#' A character vector of the column names of the baseline covariates for checking balance.

#'
#'
#' @return A \code{data.frame} containing the covariate balance measures (smd, standardized mean difference) between each cluster in the treatment arm and the control arm, both before and after the weighting adjustment.
#'
#' @export
#'
balance <- function(
    data_in,
    atekCl_results,
    covariate_names = "X_dat.1",
    ttname, Kname
){
  set.seed(124)

  cv_components <- atekCl_results$cv_components

  n <- nrow(data_in)
  J <- length(unique(data_in$K[data_in$tt==1]))
  pk.multi <- cv_components[, glue("pk.multi_pred_valid.V{1:J}")]
  pt <- cv_components$t1.x_pred_valid

  tt <- data_in[[ttname]]
  wt1 <- 1*(tt==1)/(pt)
  wt0 <- 1*(tt==0)/(1-pt)

  K <- data_in[[Kname]]
  Yfake <- data_in[, covariate_name, drop=FALSE]

  yfake.kt1.w <- map(1:J, function(k=1) {
    # y.kt1_pred <- 0
    y.kt1_w <- wt1*1*(tt==1)*(K==k)*(Yfake - 0) #+
      # y.kt1_pred * pk.dr[[k]]

    y.kt1_w
  })

  yfake.kt0.w <- map(1:J, function(k=1) {
    y.t0_w <- wt0*1*(tt==0)*pk.multi[[k]]*(Yfake-0)
    y.t0_w
  })

  diff.yfake.k <- map(1:J, function(k=1) {

    if (length(covariate_name)==1) {
      diffk <- mean(yfake.kt1.w[[k]] - yfake.kt0.w[[k]]) / sd(Yfake)

      unadj_diffk <- (mean(Yfake[K==k&tt==1])-mean(Yfake[tt==0])) / sd(Yfake)

    } else {
      diffk <- ( colMeans(yfake.kt1.w[[k]]) - colMeans(yfake.kt0.w[[k]]) ) / apply(Yfake, 2, sd)
      unadj_diffk <- (colMeans(Yfake[K==k&tt==1, ])-colMeans(Yfake[tt==0, ])) / apply(Yfake, 2, sd)
    }

    data.frame(cluster=k, covariate = names(diffk), smd=diffk, unadj_smd =unadj_diffk, row.names = NULL)
  })

  diffK <- (diff.yfake.k) %>% reduce(bind_rows)
  diffK
}
