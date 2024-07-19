#' @export

atekCl <- function(cv_folds = 4L,
                   data_in,
                   ttname = "tt", Kname = "K", Yname = "Y",
                   Xnames,
                   Yfamily = "gaussian",
                   learners_tt,
                   learners_k,
                   learners_y,
                   sensitivity = NULL,
                   Fit = "mlr"
                   ) {

  set.seed(12345)
  crossfit_res <- cluster.specific.ate(
    cv_folds = cv_folds,
    data_in = data_in,
    ttname = ttname, Kname = Kname, Yname = Yname,
    Xnames = Xnames,
    Fit = Fit, # Fit = "mlr"
    omit.tt = FALSE,
    omit.k = FALSE,
    y1model_lme = "y1k",
    Yfamily = Yfamily,
    learners_tt = learners_tt,
    learners_k = learners_k,
    learners_y = learners_y,
    combination = NULL,
    sensitivity = sensitivity
  )

  # if (nboot>0) {
  #   cv_components <- crossfit_res$cv_components
  #
  #   set.seed(12333)
  #   seeds_boot <- sample(1:1e4, nboot, replace = F)
  #   boots <- sapply(seeds_boot, boot.atek.cl,
  #                   data_in = data_in,
  #                   cv_components = cv_components,
  #                   Xnames = Xnames,
  #                   ttname = ttname,
  #                   Kname = Kname,
  #                   Yname = Yname,
  #                   Yfamily = Yfamily)
  #   crossfit_res$ate_K$ate_k_var_bootcl <- apply(boots, 1, var, na.rm=TRUE)
  #   crossfit_res$ate_K$ate_k_ci1_bootcl <- apply(boots, 1, quantile, 0.025, na.rm=TRUE)
  #   crossfit_res$ate_K$ate_k_ci2_bootcl <- apply(boots, 1, quantile, 0.975, na.rm=TRUE)
  # }

  return(crossfit_res)

}
