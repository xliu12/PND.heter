
#' @import stats
#' @import tidyverse
#' @import glue
#' @import mvtnorm
#' @import SuperLearner
#' @import nnet
#' @import ranger
#' @importFrom purrr map reduce
#' @importFrom dplyr mutate select rename all_of
#' @importFrom magrittr %>%

interact <- function( covariates, order = 2 ) {
  if(!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }
  out <- as.data.frame(model.matrix(~-1 + .^2, covariates))
  colnames(out) <- gsub(":", "_", colnames(out))
  return(out)
}

fitting.tt <- function(train_data, valid_data,
                       tmodel = "tt ~ X",
                       randomized.tt.prop = 0.5,
                       Xnames,
                       SL_library = c("glm"),
                       seed = NULL
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

    if(!is.null(seed)){
      set.seed(seed)
    }

    sl_fit <- SuperLearner(
      Y = train_data$tt,
      X = train_data[, cov_names],
      newX = valid_data[, cov_names],
      family = "binomial",
      SL.library = SL_library
    )
    # sl_pred_valid <- sl_fit$SL.predict
    sl_pred <- predict(sl_fit, newdata = valid_data[, cov_names], onlySL = F)
    if ( sum(sl_fit$coef)==0 ) {
      sl_pred$pred <- sl_pred$library.predict[,1]
    }
    sl_pred_valid <- sl_pred$pred

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
                      SL_library = c("glm"),
                      seed = NULL
) {

  # Y(1) -----



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

  # E(Y | K=k, tt=1, X)

  if( length(grep("K", ymodel)) > 0 ) {
    cov_names <- c("K", Xnames)
    y_train_data$K <- factor(y_train_data$K)

  }
  if( length(grep("K", ymodel)) ==0 ) {
    cov_names <- c(Xnames)
  }




  # with K dummies
  if (ymodel %in% c("Y ~ K + t1 + X")
  ) { # e.g., mlearners to capture interactions of K, X, W

    yformula <- formula( paste0(" ~ -1 + K + ", paste0("", c(Xnames), collapse = " + ") ) )

    y_train_X <- data.frame(model.matrix(yformula, data = y_train_data))[, -1]

    sl_fit <- SuperLearner(
      Y = y_train_data$Y,
      X = y_train_X,
      newX = y_train_X,
      family = Yfamily,
      SL.library = SL_library
    )


    sl_pred_valid_K <- lapply(unique(valid_data$K[valid_data$tt==1]), FUN = function(k=1) {
      yk_valid_X <- data.frame(model.matrix(yformula, data = y_valid_data_K[[k]]))[, -1]

      sl_pred <- predict(sl_fit, newdata = yk_valid_X, onlySL = F)
      if ( sum(sl_fit$coef)==0 ) {
        sl_pred$pred <- sl_pred$library.predict[,1]
      }
      sl_pred_valid <- sl_pred$pred


      # if( length(grep("Y.within", ymodel)) > 0 ) {
      #   y.clustermean_t1k <- cluster_means[cluster_means$K==k, "Y"]
      #   sl_pred_valid <- y.clustermean_t1k + sl_pred_valid
      # }

      return(sl_pred_valid)
    })
    sl_pred_valid <- sl_pred_valid_K

  }



  # Y(0) -----

  if (ymodel == "Y ~ t0 + X"){ # control
    cov_names <- c(Xnames)
    y_train_data <- train_data[train_data$tt==0, ]
    y_valid_data <- valid_data

    # set.seed(9999)
    if(!is.null(seed)){
      set.seed(seed)
    }

    SL_library <- SL_library
    sl_fit <- SuperLearner(
      Y = y_train_data$Y,
      X = y_train_data[, cov_names],
      newX = y_valid_data[, cov_names],
      family = Yfamily,
      SL.library = SL_library
    )
    sl_pred <- predict(sl_fit, newdata = y_valid_data[, cov_names], onlySL = F)

    if ( sum(sl_fit$coef)==0 ) {
      sl_pred$pred <- sl_pred$library.predict[,1]
    }
    sl_pred_valid <- sl_pred$pred
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
                      SL_libary = c("SL.multinom"),
                      seed = NULL
) {

  if (kmodel == "K ~ t1") {
    k_train_data <- train_data[train_data$tt==1, ]
    k_train_data$intercept <- 1
    valid_data$intercept <- 1
    cov_names <- "intercept"

    multi_fit <- SL.multinom(
      Y = (k_train_data$K - min(k_train_data$K)),
      X = k_train_data[, cov_names, drop=FALSE],
      newX = valid_data[, cov_names, drop=FALSE],
      family = "multinomial"
    )
    multi_pred_valid <- multi_fit$pred

    sl_pred_valid <- list(multi_pred_valid = multi_pred_valid)
    sl_fit <- list(multi_fit = multi_fit)

  }

  if (kmodel == "K ~ t1 + X"){
    cov_names <- c(Xnames)

    # set.seed(9999)
    if(!is.null(seed)){
      set.seed(seed)
    }
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


    if( "SL.xgboost.modified" %in% SL_libary ) {
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


    if( "SL.ranger.modified" %in% SL_libary ) {
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

    if("SL.nnet.modified" %in% SL_libary) {
      fit <- SL.nnet.modified(
        Y = as.factor(k_train_data$K),
        X = k_train_data[, cov_names],
        newX = valid_data[, cov_names]
      )
      pred_valid <- fit$pred

      sl_pred_valid <- list(multi_pred_valid = pred_valid)
      sl_fit <- list(multi_fit = fit)
    }

  }



  out <- list(
    k_pred_valid = sl_pred_valid,
    k_fit = sl_fit
  )

  return(out)
}



SL.multinom <- function (Y, X, newX, family = "multinomial", obsWeights, size = 2, ...) {

  # if( !is.factor(Y) ) {
  #   Y <- factor(Y)
  # }
  if( !is.data.frame(X) ) {
    X <- as.data.frame(X)
    if( is.null(newX) ) {
      newX <- X
    }
    newX <- as.data.frame(newX)
  }
  fit.multinom <- multinom(factor(Y) ~ ., data = X
                           # , weights = obsWeights
  )

  pred <- predict(fit.multinom, newdata = newX, type = "prob")
  # pred <- c(t(pred))

  fit <- list(object = fit.multinom)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.multinom")
  return(out)
}

SL.xgboost.modified <- function (Y, X, newX, family = "multinomial", #obsWeights,
                                 ntrees = 1000,
                                 max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                                 nthread = 1, verbose = 0, save_period = NULL, ...) {


  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }

  xgmat = xgboost::xgb.DMatrix(data = X, label = Y
                               # , weight = obsWeights
  )

  model = xgboost::xgboost(data = xgmat,
                           num_class = length(unique(Y)),
                           objective = "multi:softprob", # "multi:softmax",
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

SL.ranger.modified <- function (Y, X, newX, family = "multinomial", #obsWeights,
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
                               # case.weights = obsWeights,
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



SL.nnet.modified <- function (Y, X, newX, family="multinomial", size = 2, ...)
{

  if( is.null(newX) ) {
    newX <- X
  }
  fit.nnet <- nnet::nnet(Y ~ ., data = data.frame(Y=factor(Y), X),
                         size = size, trace = FALSE,
                         maxit = 500, linout = FALSE)
  pred <- predict(fit.nnet, newdata = newX, type = "raw")
  fit <- list(object = fit.nnet)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.nnet")
  return(out)
}



