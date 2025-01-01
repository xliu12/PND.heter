
#' Estimation of the cluster-specific treatment effects in the partially nested design.
#' @param data_in A \code{data.frame} containing all necessary variables.
#' @param ttname [\code{character}]\cr
#'  A character string of the column name of the treatment variable. The treatment variable should be dummy-coded, with 1 for the (clustered) treatment arm and 0 for the (non-clustered) control arm.
#' @param Kname [\code{character}]\cr
#'  A character string of the column name of the cluster assignment variable. This variable should be coded as 0 for individuals in the control arm, the arm without the cluster assignment.
#' @param Ynames [\code{character}]\cr
#' A character string of the column name of the outcome variable
#' @param Xnames  [\code{character}]\cr
#' A character vector of the column names of the baseline covariates.
#' @param learners_tt  [\code{character}]\cr
#' A character vector of methods for estimating the treatment model, chosen from the \code{SuperLearner} R package. Default is \code{"SL.glm"}, a generalized linear model for the binary treatment variable. Other available methods can be found using the R function \code{SuperLearner::listWrappers()}.
#' @param learners_k  [\code{character}]\cr
#' A character string of a method for estimating the cluster assignment model, which can be one of \code{"SL.multinom"} (default),  \code{"SL.xgboost.modified"}, \code{"SL.ranger.modified"}, and \code{"SL.nnet.modified"}.
#' Default is  \code{"SL.multinom"}, the multinomial regression (\code{nnet::multinom}) for the categorical cluster assignment using the treatment arm data. The other options are \code{"SL.xgboost.modified"} (gradient boosted model, \code{xgboost::xgboost}), \code{"SL.ranger.modified"} (random forest model, \code{ranger::ranger}), and \code{"SL.nnet.modified"} (neural network model, \code{"SL.nnet.modified"})  modified for fitting categorical response variable of type  multinomial.
#' @param learners_y  [\code{character}]\cr
#' A character vector of methods for estimating the outcome model, chosen from the \code{SuperLearner} R package. Default is \code{"SL.glm"}, a generalized linear model for the outcome variable, with \code{family} specified by \code{Yfamily}. Other available methods can be found using the R function \code{SuperLearner::listWrappers()}.
#' @param cv_folds [\code{numeric(1)}]\cr The number of cross-fitting folds. Default is 4.
#' @param sensitivity Specification for sensitivity parameter values on the standardized mean difference scale, which can be \code{NULL} (default) or \code{"small_to_medium"}. If \code{NULL}, no sensitivity analysis will be run. If \code{"small_to_medium"}, the function will run a sensitivity analysis for the cluster assignment ignorability assumption, and the sensitivity parameter values indicate a deviation from this assumption of magnitude 0.1 and 0.3 standardized mean difference.
#'
#'
#' @return A \code{list} containing the following components:
#'
#' \item{ate_K}{A \code{data.frame} of the estimation results.
#'
#' The columns "ate_k", "std_error", "boot_ci1", and "boot_ci2" contain the estimate, standard error estimate, and lower and upper bounds of the 0.95 confidence interval of the cluster-specific treatment effect for the cluster (indicated by column "cluster") in the same row.}
#'
#' \item{cv_components}{A \code{data.frame} of nuisance model estimates.}

#' \item{sens_results}{\code{NULL} if the argument \code{sensitivity = NULL}.
#'
#' If the argument \code{sensitivity = "small_to_medium"} is specified, \code{sens_results} is a list of four data frames, containing the estimation results with the sensitivity parameter value (standardized mean difference) being 0.1, 0.3, -0.1, -0.3.}
#'
#' @export
#'
#' @examples
#'
#' library(tidyverse)
#' library(SuperLearner)
#' library(glue)
#'
#' # data
#' data(data_in)
#' data_in <- data_in
#'
#' # baseline covariates
#' Xnames <- c(grep("X_dat", colnames(data_in), value = TRUE))
#'
#' estimates_ate_K <- atekCl(
#' data_in = data_in,
#' ttname = "tt",  # treatment variable
#' Kname = "K",    # cluster assignment variable, coded as 0 for individuals in the (non-clustered) control arm
#' Yname = "Y",    # outcome variable
#' Xnames = Xnames
#' )
#' estimates_ate_K$ate_K

#'
#'




# Estimating cluster-specific treatment effects (ate_K)--------------------------------------



atekCl <- function(data_in,
                   ttname,
                   Kname,
                   Yname,
                   Xnames,
                   Yfamily = "gaussian",
                   learners_tt = c("SL.glm"),
                   learners_k = c("SL.multinom"),
                   learners_y = c("SL.glm"),
                   sensitivity = NULL,
                   Fit = "mlr", cv_folds = 4L
) {

  data_in1 <- data_in %>%
    rename(tt = all_of(ttname), K = all_of(Kname))
  data_in1$K[data_in1$tt==0] <- 0
  data_in1$K[data_in1$tt==1] <- match(data_in1$K[data_in1$tt==1], unique(data_in1$K[data_in1$tt==1]))

  set.seed(12345)
  crossfit_res <- cluster.specific.ate(
    cv_folds = cv_folds,
    data_in = data_in1,
    ttname = "tt", Kname = "K", Yname = Yname,
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

  crossfit_res$ate_K <- crossfit_res$ate_K %>%
    mutate( cluster = unique(data_in[[Kname]][data_in[[ttname]]==1]),
            std_error = sqrt(indboot_var)) %>%
    select(cluster, ate_k, std_error, boot_ci1, boot_ci2)


  return(crossfit_res)

}


eif.k <- function(v = 1,  #fold,
                  folds,
                  data_in,
                  Xnames,
                  Fit = "lme", # Fit = "mlr"
                  omit.tt = FALSE,
                  omit.k = FALSE,
                  y1model_lme,
                  Yfamily = "gaussian",
                  cv = TRUE,
                  learners_tt,
                  learners_k,
                  learners_y
) {
  # try with one of the folds
  # fold <- origami::make_folds(data_in, fold_fun = folds_vfold, V = 2)
  # # fold[[v]]$training_set
  # train_data <- data_in[fold[[1]]$training_set, ]
  # valid_data <- data_in[fold[[1]]$validation_set, ]
  if(cv==TRUE) {
    # # make training and validation data
    # train_data <- origami::training(data_in)
    # valid_data <- origami::validation(data_in)
    train_data <- data_in[folds[[v]]$training_set, ]
    valid_data <- data_in[folds[[v]]$validation_set, ]
    # fold.ind <- origami::fold_index()
    fold.ind <- v
    valid_set <- folds[[v]]$validation_set
    valid_set_id <- data_in[valid_set, "id"]
  }
  if(cv==FALSE) {
    train_data <- valid_data <- data_in
    fold.ind <- 0
    valid_set <- 1:nrow(valid_data)
    valid_set_id <- data_in[valid_set, "id"]
  }

  # Fit models  -----------------

  tmodel <- ifelse(omit.tt == TRUE, "tt ~ 1", "tt ~ X")
  kmodel <- ifelse(omit.k == TRUE, "K ~ t1", "K ~ t1 + X")

  y1model <- ifelse(y1model_lme == "y1lme_x",
                    "Y ~ (1 | K) + t1 + X.between + X.within",
                    ifelse(y1model_lme == "y1lme_omitx", "Y ~ (1 | K) + t1",
                           "Y ~ K + t1 + X")
  )
  y1model <- ifelse(y1model_lme=="y1k_omitx", "Y ~ K + t1", y1model)

  y0model <- ifelse(y1model_lme == "y1k_omitx", "Y ~ t0",
                    "Y ~ t0 + X")

  # if parametric (generalized) linear models ----
  if (Fit == "lme") {
    t_out.x <- fitting.tt(
      train_data = train_data, valid_data = valid_data,
      tmodel = tmodel,
      SL_library = c("SL.glm"),
      Xnames = Xnames
    )

    y_out.kt1xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      SL_library = c("SL.glm"),
      ymodel = y1model,
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    y_out.t0xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      ymodel = y0model,
      SL_library = c("SL.glm"),
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    k_out.t1xw <- fitting.K(
      train_data = train_data, valid_data = valid_data,
      kmodel = kmodel,
      SL_libary = c("SL.multinom"),
      Xnames = Xnames
    )
  }

  # if use flexible regressions algorithms ----

  if (Fit != "lme") {
    t_out.x <- fitting.tt(
      train_data = train_data, valid_data = valid_data,
      tmodel = tmodel,
      SL_library = learners_tt,
      Xnames = Xnames
    )

    y_out.kt1xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      SL_library = learners_y,
      ymodel = y1model, #"Y ~ K + t1 + X",
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    y_out.t0xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      ymodel = y0model, # "Y ~ t0 + X",
      SL_library = learners_y,
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    k_out.t1xw <- fitting.K(
      train_data = train_data, valid_data = valid_data,
      kmodel = kmodel,
      SL_libary = learners_k,
      Xnames = Xnames
    )

  }

  # Estimators / eifs ----------------
  # (Jiang et al., 2022)

  # for K, principal score p(K1 = k | X) ----------------------
  # doubly robust estimator for p(k|t=1,x,w)
  t1.x_pred_valid <- bound_propensity(t_out.x$t1_pred_valid)
  k_out.t1xw$k_pred_valid$multi_pred_valid <- bound_precision(k_out.t1xw$k_pred_valid$multi_pred_valid)

  wt.x_valid <- with(valid_data, {
    as.numeric(tt==1) / ((t1.x_pred_valid))
  })
  wt0.x_valid <- with(valid_data, {
    as.numeric(tt==0) / ((1 - t1.x_pred_valid))
  })

  pk_dr_valid <- do.call(cbind,
                         lapply(unique(valid_data$K[valid_data$tt==1]),
                                FUN = function(k=1) {
                                  pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
                                  unstab_eif_pk <- with(valid_data, {
                                    wt.x_valid * as.numeric(tt==1) * (as.numeric(K==k) - pk.t1xw_pred_valid)
                                  })
                                  eif_pk <- with(valid_data, {
                                    wt.x_valid * as.numeric(tt==1) * (as.numeric(K==k) - pk.t1xw_pred_valid) / mean(wt.x_valid * as.numeric(tt==1))
                                  })
                                  dr_pk <- unstab_eif_pk + pk.t1xw_pred_valid
                                  dr_pk
                                }))

  pk_dr_valid <- as.data.frame(pk_dr_valid)


  # sensitivity -------
  pkpk_dr_valid <- do.call(cbind,
                           lapply(unique(valid_data$K[valid_data$tt==1]),
                                  FUN = function(k=1) {
                                    pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
                                    unstab_eif_pkpk <- with(valid_data, {
                                      2*pk.t1xw_pred_valid* wt.x_valid * as.numeric(tt==1) * (as.numeric(K==k) - pk.t1xw_pred_valid)
                                    })

                                    dr_pkpk <- unstab_eif_pkpk + pk.t1xw_pred_valid^2
                                    dr_pkpk
                                  }))

  pkpk_dr_valid <- as.data.frame(pkpk_dr_valid)


  # for Y, outcome mean Y(t)|K(1)=k --------------

  # y treat, cluster k
  # if (full_data != nrow(data_in)) {
  #   full_data <- rbind(train_data, valid_data)
  # }
  # if (nrow(train_data)==nrow(data_in)) {
  #   full_data <- data_in
  # }
  # cov_names <- c(Xnames)
  # cluster_means <- aggregate(full_data[full_data$tt==1, c("Y", Xnames)],
  #                            by = list(K = full_data$K[full_data$tt==1]),
  #                            mean)

  # dr: tp + ps + om -------

  # y treat
  ymean.kt1_dr_valid <- do.call(cbind,
                                lapply(unique(valid_data$K[valid_data$tt==1]),
                                       FUN = function(k=1) {

                                         y.kt1xw_pred_valid <- y_out.kt1xw$y_pred_valid[[k]]
                                         # eif for yfit
                                         unstab_eif_y.kt1 <- with(valid_data, {
                                           wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (Y - y.kt1xw_pred_valid)
                                         })
                                         eif_y.kt1 <- with(valid_data, {
                                           wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (Y - y.kt1xw_pred_valid) / mean(wt.x_valid * as.numeric(tt==1)* as.numeric(K==k) )
                                         })

                                         ymean.kt1_dr <- unstab_eif_y.kt1 + y.kt1xw_pred_valid * pk_dr_valid[, k]

                                         return(ymean.kt1_dr)

                                       }))
  ymean.kt1_dr_valid <- as.data.frame(ymean.kt1_dr_valid)

  # y control
  # valid_data_t0 <- valid_data
  # valid_data_t0$tt <- 0
  # y.t0xw_pred_valid <- predict(y_out.t0xw$y_fit, valid_data_t0[, y_out.t0xw$y_fit$varNames])$pred
  y.t0xw_pred_valid <- y_out.t0xw$y_pred_valid

  ymean.t0_dr_valid <- do.call(cbind, lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1){
    pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]

    unstab_eif_y.t0 <- with(valid_data, {
      # wt0.x_valid * as.numeric(tt==0) * pk_dr_valid[, k] * (Y - y.t0xw_pred_valid)
      wt0.x_valid * as.numeric(tt==0) * pk.t1xw_pred_valid * (Y - y.t0xw_pred_valid)
    })

    eif_y.t0 <- with(valid_data, {
      # wt0.x_valid * as.numeric(tt==0) * pk_dr_valid[, k] * (Y - y.t0xw_pred_valid) / mean(wt0.x_valid * as.numeric(tt==0) * pk_dr_valid[, k])
      wt0.x_valid * as.numeric(tt==0) * pk_dr_valid[, k] * (Y - y.t0xw_pred_valid) / mean(wt0.x_valid * as.numeric(tt==0))
    })

    ymean.t0_dr <- unstab_eif_y.t0 + y.t0xw_pred_valid * pk_dr_valid[, k]
    ymean.t0_dr
  }))
  ymean.t0_dr_valid <- as.data.frame(ymean.t0_dr_valid)



  # nuisance functions  ----------------
  y.kt1_pred_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {
             y.kt1xw_pred_valid <- y_out.kt1xw$y_pred_valid[[k]]
             return(y.kt1xw_pred_valid)

           }))
  y.kt1_pred_valid <- as.data.frame(y.kt1_pred_valid)

  pk.multi_pred_valid <- as.data.frame(k_out.t1xw$k_pred_valid$multi_pred_valid)

  colnames(pk.multi_pred_valid) <- paste0("V", 1:ncol(pk.multi_pred_valid))

  eif_out <- list(
    # fitted
    # r_out.c = r_out.c, r_out.mc = r_out.mc,
    # m_out.rc = m_out.rc, m_out.rxc = m_out.rxc,
    # y_out.mrxc = y_out.mrxc,
    # u_out.mrc = u_out.mrc,
    # v_out.rc = v_out.rc,
    # y_out.rc = y_out.rc,
    # tilting yfit
    tmle_components = data.frame(
      fold.ind = rep(v, nrow(valid_data)),
      valid_set = valid_set,
      valid_set_id = valid_set_id,
      t1.x_pred_valid = t_out.x$t1_pred_valid,
      # y dr
      ymean.kt1_dr_valid = ymean.kt1_dr_valid,
      ymean.t0_dr_valid = ymean.t0_dr_valid,
      # k dr
      pk_dr_valid = pk_dr_valid,
      # y fit
      y.kt1_pred_valid = y.kt1_pred_valid,
      y.t0xw_pred_valid = y.t0xw_pred_valid,
      # k fit
      pk.multi_pred_valid = pk.multi_pred_valid
      # sensitivity
      , pkpk_dr_valid = pkpk_dr_valid
    ),
    # fold IDs
    # fold = fold.ind
    fold = v
    # fold = origami::fold_index()
  )

  return(eif_out)
}


# Crossfit -----------------------------------------
cluster.specific.ate <- function(
    cv_folds = 2L,
    data_in,
    ttname = "tt", Kname = "K", Yname = "Y",
    Xnames,
    Fit = "lme", # Fit = "mlr"
    omit.tt = FALSE,
    omit.k = FALSE,
    y1model_lme,
    Yfamily = "gaussian",
    learners_tt,
    learners_k,
    learners_y,
    combination = NULL,
    sensitivity = NULL
) {

  data_in$id <- 1:nrow(data_in)
  data_in[[Kname]] <- as.numeric(data_in[[Kname]])
  data_in[[Kname]][data_in[[ttname]]==0] <- 0
  data_in[[Kname]][data_in[[ttname]]==1] <- match(data_in[[Kname]][data_in[[ttname]]==1], unique(data_in[[Kname]][data_in[[ttname]]==1]))
  data_in <- data_in %>%
    rename(tt = all_of(ttname), K = all_of(Kname), Y = all_of(Yname))

  if(cv_folds == 0) {
    eif_out <- eif.k(fold = 0,
                     data_in = data_in,
                     Xnames = Xnames,
                     Fit = Fit, # Fit = "mlr"
                     omit.tt = omit.tt,
                     omit.k = omit.k,
                     y1model_lme = y1model_lme,
                     Yfamily = Yfamily,
                     cv = FALSE,
                     learners_tt = learners_tt,
                     learners_k = learners_k,
                     learners_y = learners_y
    )
    cv_components <- eif_out$tmle_components
  }

  if(cv_folds > 1) {

    fold_K <- lapply(unique(data_in$K), FUN = function(k=1) {

      if (nrow(data_in[data_in$K==k, ]) >= 1) {
        fk <- origami::make_folds(data_in[data_in$K==k, ],
                                  fold_fun = origami::folds_vfold,
                                  V = cv_folds)
        fold_k <- fk
        for(v in 1:cv_folds) {
          fold_k[[v]]$validation_set <- data_in$id[data_in$K==k][fk[[v]]$validation_set]
          fold_k[[v]]$training_set <- data_in$id[data_in$K==k][fk[[v]]$training_set]
        }
      }

      # if (nrow(data_in[data_in$K==k, ]) < 4) {
      #   # if cluster size too small, no cluster split; use entire cluster as both training and valid
      #   fk <- origami::make_folds(
      #     data_in[data_in$K==k, ][sample(1:nrow(data_in[data_in$K==k, ]), cv_folds*2, replace = T), ],
      #                             fold_fun = origami::folds_vfold,
      #                             V = cv_folds
      #   )
      #   fold_k <- fk
      #   for(v in 1:cv_folds) {
      #     fold_k[[v]]$validation_set <- data_in$id[data_in$K==k]
      #     fold_k[[v]]$training_set <- data_in$id[data_in$K==k]
      #   }
      #
      # }

      return(fold_k)
    } )

    folds <- origami::make_folds(data_in,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds)

    for(v in 1:cv_folds) {
      folds[[v]]$validation_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$validation_set
      }))
      folds[[v]]$training_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$training_set
      }))
    }

    eif <- NULL
    for(v in 1:cv_folds) {
      eif_out <- eif.k(
        v = v,
        folds = folds,
        data_in = data_in,
        Xnames = Xnames,
        Fit = Fit,
        omit.tt = omit.tt,
        omit.k = omit.k,
        y1model_lme = y1model_lme,
        Yfamily = Yfamily,
        cv = TRUE,
        learners_tt = learners_tt,
        learners_k = learners_k,
        learners_y = learners_y
      )
      eif <- rbind(eif, eif_out$tmle_components)
    }
    cv_components <- eif[order(eif$valid_set), ] # make sure it is in the same order of data_in
    # sum(cv_components$valid_set - data_in$id)
  }


  # Effects -----

  ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {

    ymean.kt1 <- cv_components[, paste0("ymean.kt1_dr_valid.V", k)]
    ymean.t0 <- cv_components[, paste0("ymean.t0_dr_valid.V", k)]
    pk <- cv_components[, paste0("pk_dr_valid.V", k)]

    # point estimate
    ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

    # variance estimate
    psi_den <- mean(pk)
    # eif_k <- (ymean.kt1 - ymean.t0 - mean(ymean.kt1 - ymean.t0)) / psi_den - (ate_k/psi_den)*(pk - mean(pk)) # equal to
    eif_k <- (ymean.kt1 - ymean.t0) / psi_den - (ate_k/psi_den)*(pk)

    ate_k_var <- var(eif_k)/length(eif_k)
    wald_ci1 <- ate_k - qnorm(0.975)*sqrt(ate_k_var)
    wald_ci2 <- ate_k + qnorm(0.975)*sqrt(ate_k_var)

    c(ate_k = ate_k, ate_k_var = ate_k_var,
      wald_ci1=wald_ci1, wald_ci2=wald_ci2)
  }))
  ate_K <- as.data.frame(ate_K)

  # sensitivity -------------
  if (length(sensitivity)==0) { sens_results <- NULL }
  if (length(sensitivity) > 0) {
    if (sensitivity == "small_to_medium") {
      smdk_list <- c(-0.3, -0.1, 0.1, 0.3)
      dk_list <- smdk_list*sd(data_in[data_in$tt==0, "Y"])
    } else {
      dk_list <- smdk_list <- sensitivity
    }
    sens_results <- list()
    for (s in 1:length(dk_list)) {
      dk <- dk_list[s]

      sens_ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {

        ymean.kt1 <- cv_components[, paste0("ymean.kt1_dr_valid.V", k)]
        ymean.t0 <- cv_components[, paste0("ymean.t0_dr_valid.V", k)]
        pk <- cv_components[, paste0("pk_dr_valid.V", k)]
        pkpk <- cv_components[, paste0("pkpk_dr_valid.V", k)]

        # point estimate
        ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)
        psi_sens_k <- mean(pkpk)/mean(pk)
        SenSate_k <- ate_k - dk + dk*psi_sens_k
        # variance estimate
        psi_den <- mean(pk)
        # eif_k <- (ymean.kt1 - ymean.t0 - mean(ymean.kt1 - ymean.t0)) / psi_den - (ate_k/psi_den)*(pk - mean(pk)) # equal to
        eif_k <- (ymean.kt1 - ymean.t0) / psi_den - (ate_k/psi_den)*(pk)
        eif_sens_k <- eif_k + dk * (pkpk/psi_den - (psi_sens_k/psi_den)*pk)

        SenSate_var_k <- var(eif_sens_k)/length(eif_sens_k)
        wald_ci1 <- SenSate_k - qnorm(0.975)*sqrt(SenSate_var_k)
        wald_ci2 <- SenSate_k + qnorm(0.975)*sqrt(SenSate_var_k)


        c(SenSate_k = SenSate_k, SenSate_var_k = SenSate_var_k,
          wald_ci1=wald_ci1, wald_ci2=wald_ci2)
      }))
      sens_ate_K <- as.data.frame(sens_ate_K)

      sens_results[[glue("{smdk_list[s]}")]] <- sens_ate_K
    }

  }



  if (length(combination)>0) {
    eif_K <- lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
      ymean.kt1 <- cv_components[, paste0("ymean.kt1_dr_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_dr_valid.V", k)]
      pk <- cv_components[, paste0("pk_dr_valid.V", k)]

      # mu.kt1 <- cv_components[, paste0("y.kt1_pred_valid.V", k)]
      # mu.t0 <- cv_components[, paste0("y.t0xw_pred_valid")]
      # pk.multi <- cv_components[, paste0("pk.multi_pred_valid.V", k)]
      # plugin_den <- mean(pk.multi)
      # plugin_num <- mean(mu.kt1 - mu.t0)
      # biaseif_k <- (ymean.kt1 - ymean.t0) / plugin_den - ((plugin_num/plugin_den)/plugin_den)*(pk)
      # onestep_k <- (plugin_num/plugin_den) + mean(biaseif_k)

      # separate one-step for numerator and denominator
      ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)
      psi_den <- mean(pk)
      eif_k <- (ymean.kt1 - ymean.t0) / psi_den - (ate_k/psi_den)*(pk)

      list(eif_k=eif_k, ate_k=ate_k)
    })

    if (combination[1]=="cluster_average") {
      num_clus <- length( unique(data_in$K[data_in$tt==1]) )
      comb <- rep(1/num_clus, num_clus)
    }
    if (combination[1]!="cluster_average") {
      comb <- combination
    }
    eif_Kmat <- do.call(cbind, lapply(1:length(eif_K), function(s=1) {
      eif_K[[s]]$eif_k
    }))
    # onestep_Kmat <- do.call(cbind, lapply(1:length(eif_K), function(s=1) {
    #   eif_K[[s]]$onestep_k
    # }))
    eif_comb <- eif_Kmat %*% comb # rowMeans(eif_K)
    # ate_comb <- drop(onestep_Kmat %*% comb) # mean(ate_K$ate_k * )
    atek_comb <- sum(ate_K$ate_k * comb)
    comb_var <- var(eif_comb)/length(eif_comb)
    comb_wald_ci1 <- atek_comb - qnorm(0.975)*sqrt(comb_var)
    comb_wald_ci2 <- atek_comb + qnorm(0.975)*sqrt(comb_var)

    ateComb <- data.frame(atek_comb, comb_wald_ci1, comb_wald_ci2)
  }



  # Interval ------------------------

  set.seed(12345)
  ind_bootout <- boot::boot(data = cv_components,
                            statistic = fun.each.boot, data_in=data_in,
                            R = 1000,
                            weights = NULL)

  ind_boot_K <- ind_bootout$t
  ind_bootci <- apply(ind_boot_K, 2, FUN = function(s) { quantile(s, probs = c(0.025, 0.975), na.rm=TRUE)} )
  ind_bootvar <- apply(ind_boot_K, 2, FUN = function(s) { var(s[abs(s)<1e2], na.rm=TRUE)} )

  # collect the intervals
  ate_K$boot_ci1 <- ind_bootci[1, 1:nrow(ate_K)]
  ate_K$boot_ci2 <- ind_bootci[2, 1:nrow(ate_K)]
  ate_K$indboot_var <- ind_bootvar[1:nrow(ate_K)]

  if (length(combination)>0) {
    ateComb$bootcomb_ci1 <- ind_bootci[1, nrow(ate_K)+1]
    ateComb$bootcomb_ci2 <- ind_bootci[2, nrow(ate_K)+1]
  }
  if (length(combination)==0) {
    ateComb <- NULL
  }


  # sensitivity interval boot -------------
  if (length(sensitivity)>0) {
    set.seed(12345)

    boot_sens_results <- list()
    for (s in 1:length(dk_list)) {
      ind_bootout <- boot::boot(data = cv_components,
                                statistic = fun.each.bootsens, dk=dk_list[s],
                                data_in=data_in,
                                R = 1000,
                                weights = NULL)

      ind_boot_K <- ind_bootout$t
      ind_bootci <- apply(ind_boot_K, 2, FUN = quantile, probs = c(0.025, 0.975), na.rm=TRUE)

      sens_results[[glue("{smdk_list[s]}")]]$boot_ci1 <- ind_bootci[1, 1:nrow(ate_K)]
      sens_results[[glue("{smdk_list[s]}")]]$boot_ci2 <- ind_bootci[2, 1:nrow(ate_K)]
    }

  }


  # output ----

  crossfit_out <- list(cv_components = cv_components,
                       # dr
                       #ateComb=ateComb,
                       ate_K = ate_K
                       ,sens_results = sens_results
  )
  return(crossfit_out)
}







bound_precision <- function(vals, tol = 1e-6) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}

bound_propensity <- function(vals, bounds = c(0.01, 0.99)) {
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}
scale_from_unit <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}








