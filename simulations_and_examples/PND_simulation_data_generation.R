library(mvtnorm)

gen.data <- function(seed = 12345,
                    J=40, # number of clusters in the treatment arm
                    sample_size=40*40, # sample size
                    gen_randomized.tt = FALSE, # (randomized treatment assignment)
                    gen_quadratic.tt = FALSE, 
                    gen_quadratic.K1 = FALSE,
                    gen_quadratic.Y = FALSE
){
  
  
  ############ Data Generation ################ 
  
  # sample size 
  
  J <- J
  m <- sample_size / J
  n <- sample_size # sample size
  
  id <- 1:n
  
  # covariates 
  
  # X: influence treatment, cluster, and outcome
  
  num_x <- 4 
  
  X_dat <- rmvnorm(n, sigma = diag(1, nrow = num_x))
  
  
  # Treatment ---- 
  xmat <- model.matrix(~ -1 + ., data.frame(X_dat=X_dat))
  
  if (gen_quadratic.tt == TRUE) {
    xmat1 <- xmat
    xmat <- apply(xmat1, 2, function(x) { (x^2 - 1)/sqrt(2) })
  }
  
  alpha.X <- sqrt(0.15)
  if (gen_randomized.tt == TRUE) { alpha.X <- 0 }
  bx <- alpha.X / sqrt(ncol(xmat)) / ((1:ncol(xmat))^(0))
  
  # treated proportion
  treated_p <- 0.5
  tt <- cut(xmat %*% bx + rnorm(n), qnorm(c(0, treated_p, 1)), include.lowest = T) 
  tt <- as.numeric(tt)-1
  
  # Potential cluster assignment ---------------------
  
  wmat <- model.matrix(~ -1 + ., data.frame(X_dat=X_dat))
  
  if (gen_quadratic.K1 == TRUE) {
    wmat1 <- wmat
    wmat <- apply(wmat1, 2, function(x) { (x^2 - 1)/sqrt(2) })
  }
  bw <- sqrt((0.4) / ncol(wmat)) / ((1:ncol(wmat))^(0))
  
  K1 <- cut(wmat %*% bw + rnorm(n), qnorm((0:J)/J), include.lowest = T)
  K1 <- as.numeric(K1)

  ICC_Y1 <- 0.2
  sigmau1 <- sqrt( ICC_Y1 ) 
  sigmae1 <- sqrt( 1-ICC_Y1 )
  
  # cluster mean of baseline covariates
  cluster_K <- aggregate(data.frame(Xb_dat=X_dat), by = list(K1=K1), mean)
  
  # each cluster's characteristic 
  fixed <- TRUE
  if (fixed == TRUE) { 
    # performance for estimating/predicting cluster-specific treatment effects
    # generated cluster random effects that were then fixed across the replicates for a given scenario. While this modification may appear atypical, such an approach has previously been used when examining the estimation of random effects.
    set.seed(12345)
    # cluster intercept (not interact with individual covariates)
    cluster_K$u1 <- (rnorm(J, 0, sigmau1))
    # cluster slope (interact with individual covariates)
    cluster_K$vx1 <- (rnorm(J, 0, sigmau1/sqrt(2))) 
    cluster_K$vw1 <- (rnorm(J, 0, sigmau1/sqrt(2))) 
    set.seed(seed = seed) # within-cluster outcome residuals are still randomly varying from sample to sample
  }
  
  # data so far
  dat <- merge(data.frame(id=id, K1=K1, tt=tt, X_dat=X_dat), 
               cluster_K, by = "K1")
  
  # Potential outcomes ----
  
  ratio_ve01 <- 1
  sigmae0 <- sqrt( sigmae1^2*ratio_ve01 )
  
  # unexplained outcome variance
  revar <- sqrt( 0.5*(sigmau1^2+sigmae1^2) + 0.5*sigmae0^2 ) 
  
  ICC_X <- 0.2
  bx_control <- sqrt( revar^2*0.26/1/0.5 ) 
  bx_within <- sqrt( revar^2*0.15/(1-ICC_X)/0.5 ) 
  bx_between<- sqrt( revar^2*0.39/ICC_X/0.5 ) 
  
  e0 <- rnorm(n,0,sigmae0)
  e1 <- rnorm(n,0,sigmae1)
  
  Trteff.size <- 0.5
  gamma0 <- 0; 

  # control
  sdpool <- sqrt(bx_control^2 + bx_control^2 + sigmae0^2)
  
  gamma1 <- Trteff.size*sdpool + gamma0
  
  # treatment outcome
  
  kmat <- cbind(
    u1 = dat$u1, 
    (dat[, grep("Xb_dat", colnames(dat),value = T)] - 0)/sqrt(1),
    (dat[, grep("X_dat", colnames(dat),value = T)] - 0)/sqrt(1)
  )
  
  if (gen_quadratic.Y == TRUE) {
    kmat <- cbind(
      u1 = dat$u1, 
      vx1 = dat$vx1 * dat$X_dat.1, 
      vw1 = dat$vw1 * dat$X_dat.2,
      XbX1 = (dat$Xb_dat.1 * dat$X_dat.1 - 0)/sqrt(1), 
      (dat[, grep("Xb_dat", colnames(dat),value = T)] - 0)/sqrt(1),
      (dat[, grep("X_dat", colnames(dat),value = T)]^2 - 1)/sqrt(2*1^2)
    )
  }
  
  bk_treat <- 0 / ((1:ncol(kmat))^(0))
  names(bk_treat) <- colnames(kmat)
  
  if (gen_quadratic.Y == TRUE) {
    bk_treat[grep("vx1", colnames(kmat), value = T)] <- 1
    bk_treat[grep("vw1", colnames(kmat), value = T)] <- 1
    bk_treat[grep("XbX1", colnames(kmat), value = T)] <- bx_between * sqrt(ICC_X/ (ICC_X + ICC_X^2))
  }
  
  bk_treat[grep("u1", colnames(kmat), value = T)] <- 1
  bk_treat[grep("Xb_", colnames(kmat), value = T)] <- (bx_between-bx_within) / sqrt(2)
  bk_treat[grep("X_", colnames(kmat), value = T)] <- bx_within / sqrt(2)
  
  # control outcome
  cmat <- cbind(
    model.matrix(~ -1 + ., dat[, grep("X_dat", colnames(dat),value = T)])
  )
  if (gen_quadratic.Y == TRUE) {
    cmat <- cbind(
      (dat[, grep("X_dat", colnames(dat),value = T)]^2 - 1)/sqrt(2)
    )
  }

  bk_control <- 0 / ((1:ncol(cmat))^(0))
  names(bk_control) <- colnames(cmat)
  bk_control[grep("X", colnames(cmat), value = T)] <- bx_control / sqrt(2) 
  
  Y1 <- gamma1 + as.matrix(kmat) %*% bk_treat + e1  
  Y0 <- gamma0 + as.matrix(cmat) %*% bk_control + e0
  
  # True Values 
  ate.sample <- mean( Y1-Y0 )
  
  sample.te_K <- rep(0, J) 
  for(k in 1:J) { # E[Y1 - Y0 | K1=k]
    sample.te_K[k] <- mean(( Y1-Y0 )[dat$K1==k])
  }
  
  # Observed Data ---------------------------- 
  dat$Y <- as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt))
  dat$K <- dat$K1 * dat$tt + 0 * (1-dat$tt)
  
  datobs <- dat[, c("Y", "K", "tt", 
                    grep("X_dat", colnames(dat),value = T))]
  
  datobs$id <- 1:nrow(datobs)
  
  
  datobs <- datobs
  
  return(datobs)
}


