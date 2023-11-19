
#' Estimation of the cluster-specific treatment effects in the partially nested design
#' @param data_in A \code{data.frame} containing the observed data.
#'    In "data_in", column "tt" is the treatment assignment ("tt" is coded as 0 for individuals in the control arm and as 1 for individuals in the treatment arm);
#'    column "K" is the cluster assignment in the treatment arm ("K" is coded as 1, 2, ..., J for each individual in the treatment arm with J being the number of clusters, and "K" is coded as 0 for individuals in the control arm);
#'    column "Y" is the outcome. The other columns are baseline covariates (X).
#' @param Xnames A character vector of the names of the columns in "data_in" that correspond to baseline covariates (X)
#' @param estimator A character vector of the names of the estimators to use for estimating the cluster-specific treatment effects (ATE_k, k = 1,...,J). The estimators currently supported include those described in Liu (2023), including (i) trt-cluster, (ii) trt-y, (iii) cluster-y, (iv) triply-robust (linear), and (v) triply-robust (dml).
#'    The estimators (i)-(iv) are implemented with the parametric models where linear terms of the baseline covariates are included.
#'    the estimator (v) triply-robust (dml) is implemented with the double machine learning procedure (Chernozhukov et al., 2018 <https://doi.org/10.1111/ectj.12097>) with two-fold cross-fitting and data-adaptive packages; specifically, for estimating the cluster assignment probability, the "xgboost" R package is used <https://cran.r-project.org/package=xgboost>; for estimating the treatment probability and the outcome mean, an ensemble of algorithms is used, including boosted trees (via the “xgboost” R package), random forest (via the "ranger" R package <https://cran.r-project.org/package=ranger>), and generalized additive model (via the "gam" R package <https://cran.r-project.org/package=xgboost>.), implemented with the super learner ensembling procedure (via the "SuperLearner" R package <https://cran.r-project.org/package=SuperLearner>).
#' @param y1model_lme Whether to use the random-effects outcome regression for the estimators (ii) trt-y, (iii) cluster-y, and (iv) triply-robust (linear), where the outcome mean is involved. If "y1model_lme = TRUE", the random-intercept outcome regression (with the cluster-mean covariates and cluster-mean centered covariates) will be used; if "y1model_lme = FALSE", the fixed-effects outcome regression will be used.
#' @param randomized.tt Whether the treatment assignment is randomized. If "randomized.tt = TRUE", the treatment probability will be a constant specified by the argument "randomized.ttprop", which is the proportion of individuals randomized to the treatment arm.
#' @param randomized.ttprop The proportion of individuals randomized to the treatment arm.
#'
#'
#' @export
#'
#' @examples
#'  data(data_in)
#'  data_in <- data_in
#'  Xnames <- c(grep("X_dat", colnames(data_in), value = TRUE))
#'
#'  # estimates_ate_K <- cluster.specific.ate(
#'  # data_in = data_in,
#'  # Xnames = Xnames,
#'  # estimator = c("trt-cluster",
#'  # "trt-y",
#'  # "cluster-y",
#'  # "triply-robust (linear)",
#'  # "triply-robust (dml)"),
#'  # y1model_lme = FALSE,
#'  # randomized.tt = FALSE, randomized.ttprop = NULL
#'  # )
#'
#'




# library(lme4)
# library(SuperLearner)
# library(ranger)
# library(xgboost)
# library(nnet)
# library(SuperLearner)
# library(origami)
# library(boot)

# Estimating cluster-specific treatment effects (ate_K)--------------------------------------
cluster.specific.ate <- function(
    data_in,
    Xnames,
    estimator = c("trt-cluster", "trt-y", "cluster-y", "triply-robust (linear)", "triply-robust (dml)"),
    y1model_lme = FALSE,
    randomized.tt = FALSE, randomized.ttprop = NULL
) {

  # parametric modeling based -------------------------------------
  if (sum("triply-robust (dml)" != estimator) > 0) {
    est_out <- ate.k.one(fold = 0,
                     data_in = data_in,
                     Xnames = Xnames,
                     Fit = "lme",
                     randomized.tt = randomized.tt,
                     randomized.tt.prop = randomized.ttprop,
                     y1model_lme = y1model_lme,
                     Yfamily = "gaussian",
                     cv = FALSE
    )
    cv_components <- est_out$est_components



    # triply-robust ----
    triply_linear.ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {

      ymean.kt1 <- cv_components[, paste0("ymean.kt1_dr_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_dr_valid.V", k)]
      pk <- cv_components[, paste0("pk_dr_valid.V", k)]

      # point estimate
      ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

      c(ate_k = ate_k)
    }))
    triply_linear.ate_K <- as.data.frame(triply_linear.ate_K)


    # Interval
    # individual bootstrap scores
    dat_components <- merge(data_in, cv_components, by.x = "id", by.y = "valid_set")

    set.seed(12345)

    fun.each.boot <- function(dat_components, boot_inds) {

      boot_ate_K <- lapply(unique(dat_components$K[dat_components$tt==1]), FUN = function(k=1) {
        # bootstrap scores for estimating ate_k
        ymean.kt1 <- dat_components[boot_inds, paste0("ymean.kt1_dr_valid.V", k)]
        ymean.t0 <- dat_components[boot_inds, paste0("ymean.t0_dr_valid.V", k)]
        pk <- dat_components[boot_inds, paste0("pk_dr_valid.V", k)]

        boot_ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

        return(c(boot_ate_k = boot_ate_k))
      })
      boot_ate_K <- data.frame(do.call(rbind, boot_ate_K))

      each_boot <- c(boot_ate_k = t(boot_ate_K$boot_ate_k))
      return(each_boot)
    }

    ind_bootout <- boot(data = dat_components,
                        statistic = fun.each.boot, R = 1000,
                        weights = NULL)

    ind_boot_K <- ind_bootout$t
    ind_bootci <- apply(ind_boot_K, 2, FUN = quantile, probs = c(0.025, 0.975), na.rm=TRUE)

    # collect the intervals
    triply_linear.ate_K$ate_k_ci1 <- ind_bootci[1, ]
    triply_linear.ate_K$ate_k_ci2 <- ind_bootci[2, ]


    # trt + cluster ----

    trtclus.ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
      ymean.kt1 <- cv_components[, paste0("ymean.kt1_trtclus_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_trtclus_valid.V", k)]

      ate_k <- mean(ymean.kt1 - ymean.t0)

      c(ate_k = ate_k )
    }))
    trtclus.ate_K <- as.data.frame(trtclus.ate_K)


    # trt + y ----
    wt.x_valid <- with(data_in, {
      as.numeric(tt==1) / ((cv_components$t1.x_pred_valid))
    })
    trty.ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
      ymean.kt1 <- cv_components[, paste0("ymean.kt1_trty_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_trty_valid.V", k)]
      ate_k <- mean(ymean.kt1 - ymean.t0)
      c(ate_k = ate_k)
    }))
    trty.ate_K <- as.data.frame(trty.ate_K)

    # cluster + y ----

    clusy.ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
      ymean.kt1 <- cv_components[, paste0("ymean.kt1_clusy_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_clusy_valid.V", k)]
      pk <- cv_components[, paste0("pk.multi_pred_valid.V", k)]

      ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

      c(ate_k = ate_k)
    }))
    clusy.ate_K <- as.data.frame(clusy.ate_K)

  }

  # dml based -------------------------
  set.seed(123)
  if (sum("triply-robust (dml)" == estimator) > 0) {
    nk_min <- min(table(data_in$K))
    if (nk_min <= 2) { print("Warning: minimum cluster size is too small.")}
    cv_folds <- 2
  }
  if (cv_folds > 1) {

    fold_K <- lapply(unique(data_in$K), FUN = function(k=1) {

      if (nrow(data_in[data_in$K==k, ]) >= 1) {
        fk <- origami::make_folds(data_in[data_in$K==k, ],
                                  fold_fun = origami::folds_vfold,
                                  V = cv_folds
        )
        fold_k <- fk
        for(v in 1:cv_folds) {
          fold_k[[v]]$validation_set <- data_in$id[data_in$K==k][fk[[v]]$validation_set]
          fold_k[[v]]$training_set <- data_in$id[data_in$K==k][fk[[v]]$training_set]
        }
      }


      return(fold_k)
    } )

    folds <- origami::make_folds(data_in,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds
    )

    for(v in 1:cv_folds) {
      folds[[v]]$validation_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$validation_set
      }))
      folds[[v]]$training_set <- unlist(lapply(1:length(fold_K), FUN = function(k=1) {
        fold_K[[k]][[v]]$training_set
      }))
    }

    est_out <- NULL
    for(v in 1:cv_folds) {
      ate.k.oneout <- ate.k.one(
        v = v,
        foldsall = folds,
        data_in = data_in,
        Xnames = Xnames,
        Fit = "mlr",
        randomized.tt = randomized.tt,
        randomized.tt.prop = randomized.ttprop,
        y1model_lme = y1model_lme,
        Yfamily = "gaussian",
        cv = TRUE
      )
      est_out <- rbind(est_out, ate.k.oneout$est_components)
    }
    cv_components <- est_out[order(est_out$valid_set), ]




    # triply-robust ----
    triply_dml.ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {

      ymean.kt1 <- cv_components[, paste0("ymean.kt1_dr_valid.V", k)]
      ymean.t0 <- cv_components[, paste0("ymean.t0_dr_valid.V", k)]
      pk <- cv_components[, paste0("pk_dr_valid.V", k)]

      # point estimate
      ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

      c(ate_k = ate_k)
    }))
    triply_dml.ate_K <- as.data.frame(triply_dml.ate_K)


    # Interval
    # individual bootstrap scores
    dat_components <- merge(data_in, cv_components, by.x = "id", by.y = "valid_set")

    set.seed(12345)

    fun.each.boot <- function(dat_components, boot_inds) {

      boot_ate_K <- lapply(unique(dat_components$K[dat_components$tt==1]), FUN = function(k=1) {
        # bootstrap scores for estimating ate_k
        ymean.kt1 <- dat_components[boot_inds, paste0("ymean.kt1_dr_valid.V", k)]
        ymean.t0 <- dat_components[boot_inds, paste0("ymean.t0_dr_valid.V", k)]
        pk <- dat_components[boot_inds, paste0("pk_dr_valid.V", k)]

        boot_ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)

        return(c(boot_ate_k = boot_ate_k))
      })
      boot_ate_K <- data.frame(do.call(rbind, boot_ate_K))

      each_boot <- c(boot_ate_k = t(boot_ate_K$boot_ate_k))
      return(each_boot)
    }

    ind_bootout <- boot(data = dat_components,
                        statistic = fun.each.boot, R = 1000,
                        weights = NULL)

    ind_boot_K <- ind_bootout$t
    ind_bootci <- apply(ind_boot_K, 2, FUN = quantile, probs = c(0.025, 0.975), na.rm=TRUE)

    # collect the intervals
    triply_dml.ate_K$ate_k_ci1 <- ind_bootci[1, ]
    triply_dml.ate_K$ate_k_ci2 <- ind_bootci[2, ]


  }



  # output ---------

  if ( !"trt-cluster" %in% estimator ) { trtclus.ate_K <- NULL }
  if ( !"trt-y" %in% estimator ) { trty.ate_K <- NULL }
  if ( !"cluster-y" %in% estimator ) { clusy.ate_K <- NULL }
  if ( !"triply-robust (linear)" %in% estimator ) { triply_linear.ate_K <- NULL }
  if ( !"triply-robust (dml)" %in% estimator ) { triply_dml.ate_K <- NULL }

  estimates_out <- list(triply_dml.ate_K = triply_dml.ate_K,
                       triply_linear.ate_K = triply_linear.ate_K,
                       trtclus.ate_K = trtclus.ate_K,
                       trty.ate_K = trty.ate_K,
                       clusy.ate_K = clusy.ate_K
  )
  estimates_out <- estimates_out[unlist(lapply(estimates_out, length) > 0)]
  return(estimates_out)
}











# Calculating the scores for each individual observation -----------------

ate.k.one <- function(v = 1,  #fold,
                      foldsall = NULL,
                      data_in,
                      Xnames,
                      Fit = "lme", # Fit = "mlr"
                      randomized.tt = FALSE,
                      randomized.tt.prop = NULL,
                      y1model_lme = FALSE,
                      Yfamily = "gaussian",
                      cv = TRUE
) {

  if(cv==TRUE) {
    train_data <- data_in[foldsall[[v]]$training_set, ]
    valid_data <- data_in[foldsall[[v]]$validation_set, ]
    fold.ind <- v
    valid_set <- foldsall[[v]]$validation_set
    valid_set_id <- data_in[valid_set, "id"]
  }
  if(cv==FALSE) {
    train_data <- valid_data <- data_in
    fold.ind <- 0
    valid_set <- 1:nrow(valid_data)
    valid_set_id <- data_in[valid_set, "id"]
  }

  # Fit models for the nuisance functions -----------------

  tmodel <- ifelse(randomized.tt == TRUE, "tt ~ 1", "tt ~ X")

  y1model_lme <- ifelse(y1model_lme == TRUE,
                        "Y ~ (1 | K) + t1 + X.between + X.within", "Y ~ K + t1 + X")


  y0model_lme <- "Y ~ t0 + X"

  # if parametric (generalized) linear models ----
  if (Fit == "lme") {
    t_out.x <- fitting.tt(
      train_data = train_data, valid_data = valid_data,
      tmodel = tmodel, # "tt ~ X",
      randomized.tt.prop = randomized.tt.prop,
      SL_library = c("SL.glm"), #
      Xnames = Xnames
    )

    y_out.kt1xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      ymodel = y1model_lme,
      Xnames = Xnames,
      SL_library = c("SL.glm"),
      Yfamily = Yfamily
    )

    y_out.t0xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      ymodel = y0model_lme,
      SL_library = c("SL.glm"),
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    k_out.t1xw <- fitting.K(
      train_data = train_data, valid_data = valid_data,
      kmodel = "K ~ t1 + X",
      SL_libary = c("SL.multinom"),
      Xnames = Xnames
    )
  }

  # if use data-adaptive (e.g., machine learning) algorithms ----

  if (Fit != "lme") {
    t_out.x <- fitting.tt(
      train_data = train_data, valid_data = valid_data,
      tmodel = tmodel, # "tt ~ X",
      randomized.tt.prop = randomized.tt.prop,
      SL_library = c("SL.xgboost", "SL.ranger", "SL.gam"),
      Xnames = Xnames
    )

    y_out.kt1xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      SL_library = c("SL.xgboost", "SL.ranger", "SL.gam"), #, "SL.glmnet"
      ymodel = "Y ~ K + t1 + X",
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    y_out.t0xw <- fitting.Y(
      train_data = train_data, valid_data = valid_data,
      ymodel = "Y ~ t0 + X",
      SL_library = c("SL.xgboost", "SL.ranger", "SL.gam"),
      Xnames = Xnames,
      Yfamily = Yfamily
    )

    k_out.t1xw <- fitting.K(
      train_data = train_data, valid_data = valid_data,
      kmodel = "K ~ t1 + X",
      SL_libary = c("SL.xgboost"),
      Xnames = Xnames
    )

  }

  # Estimators  ----------------

  # doubly robust estimator for p(k|t=1,x,w)
  t1.x_pred_valid <- bound_propensity(t_out.x$t1_pred_valid)
  wt.x_valid <- with(valid_data, {
    as.numeric(tt==1) / ((t1.x_pred_valid))
  })
  wt0.x_valid <- with(valid_data, {
    as.numeric(tt==0) / ((1 - t1.x_pred_valid))
  })

  pk_dr_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {
             pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
             unstab_eif_pk <- with(valid_data, {
               wt.x_valid * as.numeric(tt==1) * (as.numeric(K==k) - pk.t1xw_pred_valid)
             })

             dr_pk <- unstab_eif_pk + pk.t1xw_pred_valid
             dr_pk
           }))

  pk_dr_valid <- as.data.frame(pk_dr_valid)


  # for Y, outcome mean Y(t)|K(1)=k --------------

  # triply robust -------

  # y treat
  ymean.kt1_dr_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {

             y.kt1xw_pred_valid <- y_out.kt1xw$y_pred_valid[[k]]
             unstab_eif_y.kt1 <- with(valid_data, {
               wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (Y - y.kt1xw_pred_valid)
             })

             ymean.kt1_dr <- unstab_eif_y.kt1 + y.kt1xw_pred_valid * pk_dr_valid[, k]

             return(ymean.kt1_dr)

           }))
  ymean.kt1_dr_valid <- as.data.frame(ymean.kt1_dr_valid)

  # y control
  valid_data_t0 <- valid_data
  valid_data_t0$tt <- 0
  y.t0xw_pred_valid <- predict(y_out.t0xw$y_fit, valid_data_t0[, y_out.t0xw$y_fit$varNames])$pred

  ymean.t0_dr_valid <- do.call(cbind, lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1){
    pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]

    unstab_eif_y.t0 <- with(valid_data, {
      wt0.x_valid * as.numeric(tt==0) * pk.t1xw_pred_valid * (Y - y.t0xw_pred_valid)
    })

    ymean.t0_dr <- unstab_eif_y.t0 + y.t0xw_pred_valid * pk_dr_valid[, k]
    ymean.t0_dr
  }))
  ymean.t0_dr_valid <- as.data.frame(ymean.t0_dr_valid)


  # trt + cluster -------

  # y treat
  ymean.kt1_trtclus_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {
             pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
             est_y.kt1 <- with(valid_data, {

               wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (Y) / mean(pk.t1xw_pred_valid)

             })

             return(est_y.kt1)

           }))
  ymean.kt1_trtclus_valid <- as.data.frame(ymean.kt1_trtclus_valid)

  # y control
  ymean.t0_trtclus_valid <- do.call(cbind, lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1){
    pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
    est_y.t0 <- with(valid_data, {
      wt0.x_valid * as.numeric(tt==0) * pk.t1xw_pred_valid * (Y) / mean(pk.t1xw_pred_valid)

    })

    return(est_y.t0)
  }))
  ymean.t0_trtclus_valid <- as.data.frame(ymean.t0_trtclus_valid)


  # trt + y --------

  # y treat
  ymean.kt1_trty_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {

             y.kt1xw_pred_valid <- y_out.kt1xw$y_pred_valid[[k]]
             est_y.kt1 <- with(valid_data, {
               wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (y.kt1xw_pred_valid) / mean(wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) )
             })

             return(est_y.kt1)

           }))
  ymean.kt1_trty_valid <- as.data.frame(ymean.kt1_trty_valid)

  # y control
  valid_data_t0 <- valid_data
  valid_data_t0$tt <- 0
  y.t0xw_pred_valid <- predict(y_out.t0xw$y_fit, valid_data_t0[, y_out.t0xw$y_fit$varNames])$pred

  ymean.t0_trty_valid <- do.call(cbind, lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1){
    est_y.t0 <- with(valid_data, {
      wt.x_valid * as.numeric(tt==1) * as.numeric(K==k) * (y.t0xw_pred_valid) / mean(wt.x_valid * as.numeric(tt==1) * as.numeric(K==k))
    })

    return(est_y.t0)
  }))
  ymean.t0_trty_valid <- as.data.frame(ymean.t0_trty_valid)



  # cluster + y ----------

  # y treat
  ymean.kt1_clusy_valid <- do.call(
    cbind,
    lapply(unique(valid_data$K[valid_data$tt==1]),
           FUN = function(k=1) {

             pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
             y.kt1xw_pred_valid <- y_out.kt1xw$y_pred_valid[[k]]
             est_y.kt1 <- with(valid_data, {
               pk.t1xw_pred_valid * y.kt1xw_pred_valid
             })

             return(est_y.kt1)

           }))
  ymean.kt1_clusy_valid <- as.data.frame(ymean.kt1_clusy_valid)

  # y control
  valid_data_t0 <- valid_data
  valid_data_t0$tt <- 0
  y.t0xw_pred_valid <- predict(y_out.t0xw$y_fit, valid_data_t0[, y_out.t0xw$y_fit$varNames])$pred

  ymean.t0_clusy_valid <- do.call(cbind, lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1){
    pk.t1xw_pred_valid <- k_out.t1xw$k_pred_valid$multi_pred_valid[, k]
    est_y.t0 <- with(valid_data, {
      pk.t1xw_pred_valid * y.t0xw_pred_valid
    })
    return(est_y.t0)
  }))
  ymean.t0_clusy_valid <- as.data.frame(ymean.t0_clusy_valid)



  # nuisance functions ----------------
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

  est_out <- list(

    est_components = data.frame(
      fold.ind = rep(v, nrow(valid_data)),
      valid_set = valid_set,
      valid_set_id = valid_set_id,
      t1.x_pred_valid = t_out.x$t1_pred_valid,
      # triply robust
      ymean.kt1_dr_valid = ymean.kt1_dr_valid,
      ymean.t0_dr_valid = ymean.t0_dr_valid,
      # trtclus
      ymean.kt1_trtclus_valid = ymean.kt1_trtclus_valid,
      ymean.t0_trtclus_valid = ymean.t0_trtclus_valid,
      # trty
      ymean.kt1_trty_valid = ymean.kt1_trty_valid,
      ymean.t0_trty_valid = ymean.t0_trty_valid,
      # clusy
      ymean.kt1_clusy_valid = ymean.kt1_clusy_valid,
      ymean.t0_clusy_valid = ymean.t0_clusy_valid,
      # k dr
      pk_dr_valid = pk_dr_valid,
      # y fit
      y.kt1_pred_valid = y.kt1_pred_valid,
      y.t0xw_pred_valid = y.t0xw_pred_valid,
      # k fit
      pk.multi_pred_valid = pk.multi_pred_valid
    ),
    # fold IDs
    fold = v
  )

  return(est_out)
}



# Fitting nuisance functions --------------------------------------------------

fitting.tt <- function(train_data, valid_data,
                       tmodel = "tt ~ X",
                       randomized.tt.prop = NULL,
                       Xnames,
                       SL_library = c("glm")
) {


  if (tmodel == "tt ~ 1"){ # randomized
    cov_names <- c(1)
    sl_fit <- NULL
    sl_pred_valid <- rep(randomized.tt.prop, nrow(valid_data))
  }


  if (tmodel != "tt ~ 1") {
    if (tmodel == "tt ~ X"){
      cov_names <- c(Xnames)
    }

    set.seed(9999)

    sl_fit <- SuperLearner(
      Y = train_data$tt,
      X = train_data[, cov_names],
      newX = valid_data[, cov_names],
      family = "binomial",
      SL.library = SL_library
    )
    sl_pred_valid <- sl_fit$SL.predict

  }

  out <- list(
    t1_pred_valid = sl_pred_valid,
    t1_fit = sl_fit
  )

  return(out)
}



fitting.Y <- function(train_data, valid_data,
                      ymodel,
                      Xnames,
                      Yfamily = "gaussian",
                      SL_library = c("glm")
) {

  # Y(1) -----

  if (sum(colMeans(train_data) - colMeans(valid_data)) == 0) { full_data <- train_data }
  if (sum(colMeans(train_data) - colMeans(valid_data)) != 0) { full_data <- rbind(train_data, valid_data) }

  cluster_means <- aggregate(full_data[full_data$tt==1, c("Y", Xnames)],
                             by = list(K = full_data$K[full_data$tt==1]),
                             mean)

  if (length(grep("within", ymodel)) == 0) {
    y_train_data <- train_data[train_data$tt==1, ]
    y_valid_data_K <- lapply(unique(valid_data$K[valid_data$tt==1]),
                             FUN = function(k=1) {
                               valid_data_t1k <- valid_data

                               valid_data_t1k$tt <- 1
                               valid_data_t1k$K <- k

                               valid_data_t1k$K <- factor(valid_data_t1k$K, levels = levels(factor(y_train_data$K)))

                               return(valid_data_t1k)
                             })
  }

  # E(Y | K=k, tt=1, X, W)
  if ( length(grep("X.within", ymodel)) > 0 & length(grep("W.within", ymodel)) > 0 ) {
    if ( length(grep("Y.within", ymodel)) > 0 ) {
      within_names <- c("Y", Xnames)
    }
    if ( length(grep("Y.within", ymodel)) == 0 ) {
      within_names <- c(Xnames)
    }
    # train data
    y_train_data <- train_data[train_data$tt==1, ]
    for(k in unique(y_train_data$K)) {
      y_train_data[y_train_data$K == k, within_names] <- sweep(
        y_train_data[y_train_data$K == k, within_names],
        2,
        STATS = as.numeric(cluster_means[cluster_means$K==k, within_names]),
        FUN = `-`)
    }

    # valid data

    y_valid_data_K <- lapply(unique(valid_data$K[valid_data$tt==1]),
                             FUN = function(k=1) {
                               valid_data_t1k <- valid_data

                               valid_data_t1k[, c(Xnames)] <- sweep(
                                 valid_data_t1k[, c(Xnames)], 2,
                                 STATS = as.numeric(cluster_means[cluster_means$K==k, c(Xnames)]),
                                 FUN = `-`)
                               valid_data_t1k$tt <- 1
                               valid_data_t1k$K <- k

                               valid_data_t1k$K <- factor(valid_data_t1k$K, levels = levels(factor(y_train_data$K)))

                               return(valid_data_t1k)
                             })

  }

  if( length(grep("K", ymodel)) > 0 ) {
    cov_names <- c("K", Xnames)
    y_train_data$K <- factor(y_train_data$K)

  }
  if( length(grep("K", ymodel)) ==0 ) {
    cov_names <- c(Xnames)
  }



  # lmer
  if (ymodel == "Y ~ (1 | K) + t1 + X.between + X.within") {
    y_train_data_wb <- merge(y_train_data, cluster_means, by = "K", suffixes = c("", "_between"))
    yformula <- formula(paste0("Y ~ (1 | K) + ", paste0(c(Xnames), "_between", collapse = " + "), "+", paste0(c(Xnames), "", collapse = " + ")))
    if(Yfamily == "gaussian") {
      sl_fit <- lmer(yformula, data = y_train_data_wb)
    }
    if(Yfamily != "gaussian") {
      sl_fit <- glmer(yformula, data = y_train_data_wb, family = Yfamily)
    }

    # predictions
    sl_pred_valid_K <- lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1) {
      yk_valid_data_wb <- merge(y_valid_data_K[[k]], cluster_means[cluster_means$K==k, ], by = "K", suffixes = c("", "_between"))

      sl_pred_valid <- predict(sl_fit, yk_valid_data_wb, type="response")

      return(sl_pred_valid)
    })
    sl_pred_valid <- sl_pred_valid_K

  }

  # with K dummies
  if (ymodel %in% c("Y ~ K + t1 + X")
  ) { # e.g., mlearners to capture interactions of K, X, W

    yformula <- formula( paste0(" ~ -1 + K + ", paste0("", c(Xnames), collapse = " + ") ) )

    y_train_X <- data.frame(model.matrix(yformula, data = y_train_data))

    suppressWarnings(
      sl_fit <- SuperLearner(
        Y = y_train_data$Y,
        X = y_train_X,
        newX = y_train_X,
        family = Yfamily,
        SL.library = SL_library
      )
    )


    sl_pred_valid_K <- lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1) {
      yk_valid_X <- data.frame(model.matrix(yformula, data = y_valid_data_K[[k]]))

      if( length(grep("Y.within", ymodel)) == 0 ) {
        sl_pred_valid <- predict(sl_fit, yk_valid_X)$pred
      }
      if( length(grep("Y.within", ymodel)) > 0 ) {
        y.clustermean_t1k <- cluster_means[cluster_means$K==k, "Y"]
        sl_pred_valid <- y.clustermean_t1k + predict(sl_fit, yk_valid_X)$pred
      }

      return(sl_pred_valid)
    })
    sl_pred_valid <- sl_pred_valid_K

  }



  # Y(0) -----

  if (ymodel == "Y ~ t0 + X"){ # control
    cov_names <- c(Xnames)
    y_train_data <- train_data[train_data$tt==0, ]
    y_valid_data <- valid_data

    set.seed(9999)

    SL_library <- SL_library
    sl_fit <- SuperLearner(
      Y = y_train_data$Y,
      X = y_train_data[, cov_names],
      newX = y_valid_data[, cov_names],
      family = Yfamily,
      SL.library = SL_library
    )
    sl_pred_valid <- sl_fit$SL.predict
  }

  out <- list(
    y_pred_valid = sl_pred_valid,
    y_fit = sl_fit
  )

  return(out)
}



fitting.K <- function(train_data, valid_data,
                      kmodel = "K ~ t1 + X",
                      Xnames,
                      obs_weights = "obs_weights_Y",
                      SL_libary = c("SL.multinom")
) {



  if (kmodel == "K ~ t1 + X"){
    cov_names <- c(Xnames)

    set.seed(9999)

    k_train_data <- train_data[train_data$tt==1, ]

    if( "SL.multinom"  %in% SL_libary ) {
      multi_fit <- SL.multinom(
        Y = (k_train_data$K - min(k_train_data$K)),
        X = k_train_data[, cov_names],
        newX = valid_data[, cov_names],
        family = "multinomial"
      )
      multi_pred_valid <- multi_fit$pred

      sl_pred_valid <- list(multi_pred_valid = multi_pred_valid)
      sl_fit <- list(multi_fit = multi_fit)
    }


    if( "SL.xgboost" %in% SL_libary ) {
      xgb_fit <- SL.xgboost.modified(
        Y = (k_train_data$K - min(k_train_data$K)),
        X = k_train_data[, cov_names],
        newX = valid_data[, cov_names],
        family = "multinomial"
      )
      xgb_pred_valid <- xgb_fit$pred

      sl_pred_valid <- list(multi_pred_valid = xgb_pred_valid)
      sl_fit <- list(multi_fit = xgb_fit)
    }


    if( "SL.ranger" %in% SL_libary ) {
      ranger_fit <- SL.ranger.modified(
        Y = (k_train_data$K - min(k_train_data$K)),
        X = k_train_data[, cov_names],
        newX = valid_data[, cov_names],
        family = "multinomial"
      )
      ranger_pred_valid <- ranger_fit$pred

      sl_pred_valid <- list(multi_pred_valid = ranger_pred_valid)
      sl_fit <- list(multi_fit = ranger_fit)
    }

  }



  out <- list(
    k_pred_valid = sl_pred_valid,
    k_fit = sl_fit
  )

  return(out)
}




SL.multinom <- function (Y, X, newX, family = "multinomial", obsWeights, size = 2, ...) {

  if( !is.data.frame(X) ) {
    X <- as.data.frame(X)
    if( is.null(newX) ) {
      newX <- X
    }
    newX <- as.data.frame(newX)
  }
  suppressMessages(
    fit.multinom <- multinom(factor(Y) ~ ., data = X, trace = F)
  )


  pred <- predict(fit.multinom, newdata = newX, type = "prob")

  fit <- list(object = fit.multinom)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.multinom")
  return(out)
}

SL.xgboost.modified <- function (Y, X, newX, family = "multinomial",
                                 ntrees = 1000,
                                 max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                                 nthread = 1, verbose = 0, save_period = NULL, ...) {
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }

  xgmat = xgboost::xgb.DMatrix(data = X, label = Y)

  model = xgboost::xgboost(data = xgmat,
                           num_class = length(unique(Y)),
                           objective = "multi:softprob",
                           nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                           eta = shrinkage, verbose = verbose,
                           nthread = nthread, params = params, save_period = save_period,
                           eval_metric = "mlogloss")

  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  pred = matrix(pred, ncol = length(unique(Y)), byrow = T)

  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)

}

SL.ranger.modified <- function (Y, X, newX, family = "multinomial",
                                num.trees = 500, mtry = floor(sqrt(ncol(X))),
                                write.forest = TRUE,
                                probability = TRUE,
                                min.node.size = 10, # default
                                replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                                num.threads = 1, verbose = T, ...) {
  if( !is.data.frame(X) ) {
    X <- as.data.frame(X)
    if( is.null(newX) ) {
      newX <- X
    }
    newX <- as.data.frame(newX)
  }
  fit_ranger <- ranger::ranger(Y ~ ., data = cbind(Y = factor(Y), X),
                               num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                               replace = replace, sample.fraction = sample.fraction,
                               write.forest = write.forest,
                               probability = TRUE,
                               num.threads = num.threads,
                               verbose = verbose)
  pred <- predict(fit_ranger, data = newX)$predictions

  fit <- list(object = fit_ranger, verbose = verbose)

  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}










bound_precision <- function(vals, tol = 1e-4) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}

bound_propensity <- function(vals, bounds = c(0.01, 0.99)) {
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}







