library(mvtnorm)
library(tidyverse)
library(glue)
library(lme4)
library(SuperLearner)
library(ranger)
library(glmnet)
library(lightgbm)
library(xgboost)
library(hal9001)
#library(e1071)
library(nnet)
library(origami)
library(boot)


GenData <- function(seed = 12345,
                    J=40, # number of Kters in treatment condition
                    m=40, # potential Kter size
                    gen_randomized.tt = FALSE, # alpha.X=0 (randomized treatment assignment)
                    gen_randomized.K1 = FALSE,
                    gen_quadratic.tt = FALSE,
                    gen_quadratic.K1 = FALSE,
                    gen_quadratic.Y = FALSE,
                    te_k_homo=FALSE, # te_k_homo=TRUE (shortcut to set homogeneous treatment effects across the clusters K1)
                    if_interference=FALSE,
                    if_size=FALSE
){

  set.seed(seed=seed+1)

  # sample size
  # J <- 10
  # m <- 30
  J <- J
  m <- m
  n <- J*m # sample size

  id <- 1:n

  # Covariates ----

  # X: influence treatment, cluster, and outcome

  num_x <- 4

  covariance_x <- diag(1, nrow = num_x)
  covariance_x[covariance_x==0] <- 0.2
  if (fix_x) {
    set.seed(12345)
    X_dat <- rmvnorm(n, sigma = covariance_x)
    set.seed(seed=seed+1)
  } else {
    X_dat <- rmvnorm(n, sigma = covariance_x)
  }
  X_dat[,4] <- 1*(X_dat[,4]>0) # dummy indicator of gender

  dat <- data.frame(id=id, X_dat=X_dat)

  # Treatment ----
  xmat <- model.matrix(~ -1 + ., data.frame(X_dat=X_dat))
  if (gen_quadratic.tt == TRUE) {
    xmat <- apply(xmat, 2, function(x) { (x^2 - 1)/sqrt(2) })
  }

  alpha.X <- sqrt(0.15)
  if (gen_randomized.tt == TRUE) { alpha.X <- 0 }
  bx <- alpha.X / sqrt(ncol(xmat)) / ((1:ncol(xmat))^(0))


  tt <- 1*(xmat %*% bx + rnorm(n) >0)

  dat$tt <- as.numeric(tt)

  # Cluster assignment ---------------------

  wmat <- model.matrix(~ -1 + ., data.frame(X_dat=X_dat))
  if (gen_quadratic.K1 == TRUE) {
    wmat <- apply(wmat, 2, function(x) { (x^2 - 1)/sqrt(2) })
  }
  # intercept to set different cluster sizes
  # different clusters have different relationships with the covariates
  bw <- matrix(0, ncol = J, nrow = ncol(wmat)+1)
  colnames(bw) <- paste0("K", 1:J)
  rownames(bw) <- c("intercept", colnames(wmat))
  bw["intercept", ] <- rep(c(0, 0.05, 0.1), length.out=J) # small to large clusters
  bw[2, ] <- seq(0, 0.15, length.out=J)
  bw[3, ] <- seq(0, -0.15, length.out=J)
  bw[4, ] <- rep(c(0, 0.1), length.out=J)
  bw[5, ] <- rep(c(0, -0.1), length.out=J)

  bw[2:5, 1] <-0
  strong.K <- TRUE
  if (strong.K) {
    set.seed(12345)
    bw[2, 2:J] <- rep(seq(0.3, 2, by=0.2), length.out=J-1)
    bw[3, 2:J] <- rep(rev(seq(0.3, 2, by=0.2)), length.out=J-1)
    bw[4, 2:J] <- rep(c(-1, 1), length.out=J-1)
    bw[5, 2:J] <- rep(c(1, -1), length.out=J-1)
    set.seed(seed = seed+1)
  }

  if (gen_randomized.K1) {
    bw[2:5, ] <- 0
  }
  bw <- bw - bw[, 1] # cluster 1 is the reference

  p <- (function(h) h / rowSums(h))(exp(cbind(1, wmat) %*% bw))
  mChoices <- t(apply(p, 1, rmultinom, n=1, size=1))

  K1 <- apply(mChoices, 1, function(x) {which(x==1)})

  dat$K1 <- K1

  dat$tt <- as.numeric(tt)

  ## cluster-specific intercepts/slopes ----
  ICC_Y1 <- 0.2
  sigmau1 <- sqrt( ICC_Y1 )
  sigmae1 <- sqrt( 1-ICC_Y1 )

  cluster_K <- data.frame(K1=1:J)
  cluster_K$capacity <- 0.5*bw["intercept", ]

  ## random_uv 
  random_uv <- TRUE # randomness of clusters, and fixed factors
  if (random_uv) {
    set.seed(123456)
    cluster_K$u1 <- 0.3*c(rep(0, each=round(J/2)), rep(1, J-round(J/2))) # teacher's gender
    cluster_K$vx1 <- runif(J, 0, 0) #  below have random slope
    set.seed(seed = seed)
    # fixed
    cluster_K$capacity <- 0.3*bw["intercept", ] + 0.3*bw[2, ] # teacher's experience also affects the selection
    # fixed + random
    cluster_K$u1 <- cluster_K$u1 + rnorm(J, 0, sigmau1)
    cluster_K$vx1 <- cluster_K$vx1 + rnorm(J, 0, sigmau1/sqrt(2))
  }




  # data so far
  dat1 <- data.frame(id=id, X_dat=X_dat, K1=K1, tt=tt) %>%
    mutate(K = K1*tt + 0*(1-tt)) %>%
    group_by(K1) %>%
    mutate(across(.cols = starts_with("X_dat"),
                  .fns = list(Xb_dat = ~1*(tt==0)*( sum(.x*(tt==1))/sum(tt==1) ) + 1*(tt==1)*( (sum(.x*(tt==1))-.x)/(sum(tt==1)-1) ) ),
                  # .fns = list(Xb_dat = ~( sum(.x*(tt==1))/sum(tt==1) ) ), # very similar results to this above
                  .names = "b{.col}")) %>%
    # mutate(obs_size = sum(tt==1) / (2*n/J) ) # very similar results to this below
    mutate(obs_size = ((tt==0)*sum(tt==1) + (sum(tt==1)-1)*(tt==1)) / (2*n/J) )

  colnames(dat1)[str_detect(colnames(dat1), "bX_dat")] <- glue("Xb_dat.{1:ncol(X_dat)}")

  if (if_interference==FALSE) {
    dat1[, c(glue("Xb_dat.{1:ncol(X_dat)}"), "obs_size")] <- 0
  }

  dat1 <- dat1 %>%
    full_join(cluster_K[, c("K1", "u1", "vx1", "capacity")], by = "K1") %>%
    ungroup() %>%
    arrange(K1)

  dat <- dat1
  
  # Outcome ----
  # Y(1)|K(1)=k

  # Y(0)
  ratio_ve01 <- 1
  sigmae0 <- sqrt( sigmae1^2*ratio_ve01 )

  # unexplained outcome variance
  revar <- 1
  # X: medium effect size f^2<-0.15; var(X)<-1
  ICC_X <- 0.2
  bx_control <-sqrt( revar^2*0.26/1/0.5 ) #
  bx_within <- sqrt( revar^2*0.15/(1-ICC_X)/0.5 )
  bx_between<- sqrt( revar^2*0.26/ICC_X/0.5 )

  e0 <- rnorm(n,0,sigmae0)
  e1 <- rnorm(n,0,sigmae1)

  # treatment effect: standardized mean difference (Cohen's d)
  # medium Cohen's d<-0.5;
  Trteff.size<-1
  gamma0 <- 0;
  gamma1 <- Trteff.size + gamma0

  # treatment outcome

  kmat <- cbind(
    dat[, c("u1", "capacity", "obs_size")],
    vx1 = dat$vx1 * dat$X_dat.1,
    XbX1 = dat$Xb_dat.1 * dat$X_dat.1,
    dat[, glue("Xb_dat.{1:ncol(X_dat)}")],
    dat[, glue("X_dat.{1:ncol(X_dat)}")]
  )

  if (gen_quadratic.Y == TRUE) {
    kmat[, glue("X_dat.{1:ncol(X_dat)}")] <- (dat[, glue("X_dat.{1:ncol(X_dat)}")]^2 - 1)/sqrt(2) # quadratic terms

  }

  bk_treat <- 0 / ((1:ncol(kmat))^(0))
  names(bk_treat) <- colnames(kmat)

  bk_treat[grep("vx1", colnames(kmat), value = T)] <- 1
  bk_treat[c("u1")] <- c(1)
  bk_treat[c("capacity")] <- 1*(if_size)
  bk_treat[glue("X_dat.{1:ncol(X_dat)}")] <- bx_within / sqrt(num_x)

  bk_treat[c("obs_size")] <- 1*(if_size)
  bk_treat[glue("Xb_dat.{2:ncol(X_dat)}")] <- (if_interference)*(bx_between-bx_within) / sqrt(num_x)
  bk_treat[c("XbX1")] <- (if_interference)*bx_within / sqrt(num_x)

  # control outcome
  cmat <- cbind(
    model.matrix(~ -1 + ., dat[, glue("X_dat.{1:ncol(X_dat)}")])
  )
  if (gen_quadratic.Y == TRUE) {
    cmat <- (dat[, glue("X_dat.{1:ncol(X_dat)}")]^2 - 1)/sqrt(2)
  }

  bk_control <- 0 / ((1:ncol(cmat))^(0))
  names(bk_control) <- colnames(cmat)
  bk_control[glue("X_dat.{1:ncol(X_dat)}")] <- bx_control / sqrt(num_x)

  # homogeneous treatment effects across clusters
  if (te_k_homo == TRUE) {
    bk_treat[grep("Xb", colnames(kmat), value = T)] <- 0
    bk_treat[grep("X_", colnames(kmat), value = T)] <- bk_control[grep("X_", colnames(cmat), value = T)]
    bk_treat[c("u1", "vx1", "capacity", "obs_size")] <- 0
  }

 
  Y1 <- gamma1 + as.matrix(kmat) %*% bk_treat + e1
  Y0 <- gamma0 + as.matrix(cmat) %*% bk_control + e0


  # True Values ---------------------------------

  sample.te_K <- rep(0, J)
  sample.te_Knum <- sample.te_Kden <- rep(0, J)

  for(k in 1:J) { # E[Y1 - Y0 | K1=k]
    sample.te_K[k] <- mean(( Y1-Y0 )[dat$K1==k])
  }

  # Observed Data ----------------------------
  # dat[1:3, ]
  dat$Y <- as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt))

  sizeKobs <- (table(dat$K[dat$tt==1]))

  datobs <- dat[, c("Y", "K", "tt",
                    grep("X_dat", colnames(dat),value = T))]

  datobs$id <- 1:nrow(datobs)


  out <- mget( ls(), envir = environment() )

  return(out)
}


# Simulation----------------------------

condition_all <- map2_df(
  .x = c(10, 10, 16, 20, 20),
  .y = c(30, 40, 25, 40, 30),
  .f = function(.x, .y) {
    data.frame(expand.grid(
      J = .x,
      m = .y,
      # generating data
      gen_quadratic.tt = c(FALSE, TRUE),
      gen_quadratic.K1 = c(FALSE, TRUE),
      gen_quadratic.Y = c(FALSE, TRUE),
      gen_randomized.tt = c(FALSE, TRUE),
      gen_randomized.K1 = c(FALSE, TRUE),
      te_k_homo = c(FALSE, TRUE),
      if_interference = c(TRUE, FALSE),
      if_size = c(TRUE, FALSE)
    ))
  }) 


condition <- condition_all %>%
  filter((
         (
             (if_interference+if_size ==2 & te_k_homo==0) |
            (if_interference+if_size ==0 & te_k_homo==1) |
                 (if_interference+if_size ==0 & te_k_homo==0)
          )
         ) %>%
  # te_k_homo = 1 implies if_interference+if_size==0
  filter(gen_quadratic.tt+ gen_quadratic.Y+ gen_quadratic.K1==3)


datseeds <- c(sample(1:1e6, 1000), sample(333:3e6, 1000))

iseed <-1
cond <- 1


OneData <- function(iseed = 1, cond = 1) {
  nk_min <- 0
  seedone <- datseeds[iseed]
  while(nk_min <= 3) { # ensure at least nk>=3 to estimate ATEk, criteria of PND (Lohr 2014)
    gen_data <- GenData(seed = seedone,
                        J = condition$J[cond],
                        m = condition$m[cond],
                        # randomized
                        gen_randomized.tt = condition$gen_randomized.tt[cond],
                        gen_randomized.K1 = condition$gen_randomized.K1[cond],
                        # quadratic
                        gen_quadratic.tt = condition$gen_quadratic.tt[cond],
                        gen_quadratic.K1 = condition$gen_quadratic.K1[cond],
                        gen_quadratic.Y = condition$gen_quadratic.Y[cond],
                        # homogeneous ate_k
                        te_k_homo = condition$te_k_homo[cond],
                        if_interference = condition$if_interference[cond],
                        if_size = condition$if_size[cond]
    )
    data_tmp <- gen_data$datobs
    nk_min <- min(table(data_tmp$K))
    seedone <- seedone + 1e6
  }

  data_in <- gen_data$datobs

  ttname <- grep("^tt", colnames(data_in), value = T)
  Kname <- grep("^K", colnames(data_in), value = T)
  Yname <- grep("^Y", colnames(data_in), value = T)
  Xnames <- grep("X_dat", colnames(data_in), value = T)

  meth <- 1
  one.method <- function (meth=1) {
    
    crossfit_res <- atekCl(data_in = data_in,
                       ttname = ttname, Kname = Kname, Yname = Yname,
                       Xnames = Xnames,
                       Yfamily = "gaussian",
                       learners_tt = c("SL.nnet", "SL.ranger"),
                       learners_k = c("SL.nnet.modified"),
                       learners_y = c("SL.nnet", "SL.ranger") )
    true_te_K <- gen_data$sample.te_K

    crossfit_te_K <- data.frame(condition[cond, ], 
                                cluster_K = unique(data_in$K[data_in$tt==1]),
                                true_te_K = true_te_K,
                                crossfit_res$ate_K,
                                sizeKobs = table(data_in$K[data_in$tt==1]),
                                row.names = NULL)


    crossfit_te_K
  }

  estimates <- purrr::map_df(.x = 1, one.method)


  out <- list(
    crossfit_te_K = estimates
   , rangeKobs = range(gen_data$sizeKobs)
  )

  return(out)
}




