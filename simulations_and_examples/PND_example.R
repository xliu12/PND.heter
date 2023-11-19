
load("covariates_data.RData")


# install or download the package for implementing the estimators
# library(devtools)
# install_github("xliu12/PND.heter/PND.heter.cluster")
library(PND.heter.cluster)
# the R function information
?cluster.specific.ate


# List the covariates
covariateNames <- c(#student-level covariates:
  "X1RTHETK1", #X1 READING IRT THETA SCORE-K1 DATA FILE
  "X1MTHETK1", #X1 MATH IRT THETA SCORE--K1 DATA FILE
  "X1TCHAPP", #X1 TEACHER REPORT APPROACHES TO LEARNING
  "X1TCHCON", #X1 TEACHER REPORT SELF-CONTROL
  "X1TCHPER", #X1 TEACHER REPORT INTERPERSONAL
  "X1TCHEXT", #X1 TEACHER REPORT EXTERN PROB BEHAVIORS
  "X1TCHINT", #X1 TEACHER REPORT INTERN PROB BEHAVIORS
  "X1ATTNFS", #X1 TEACHER REPORT ATTENTIONAL FOCUS
  "X1INBCNT", #X1 TEACHER REPORT INHIBITORY CONTROL
  "X_CHSEX_R",#Sex
  "X_RACETH_R",#Race
  "X2DISABL2", #Disability
  "C1ENGHM", #English Language Learner
  "X12MOMAR", #Parents married
  "X1NUMSIB", # of siblings
  "P1OLDMOM", #age of mom
  "P1CHLDBK", #books at home
  "P2DISTHM", #distance from school
  "T2PARIN", #parental involvement
  "X12PAR1ED_I", #mother education
  "X12PAR2ED_I", #father education
  "X2INCCAT_I", #income
  "X1PAR1EMP" #PARENT 1 EMPLOYMENT STATUS
)


# simulate treatment assignment, cluster assignment, outcome given the covariates
# using data-generation models similar to Scenario 1 of the simulation study

# packages may be needed
library(mvtnorm)
library(SuperLearner)
library(lme4)
library(xgboost)
library(ranger)
library(nnet)
library(origami)
library(boot)

real.mimic <- function(J=34, n = (470+823), te_k_homo = F){

  n_treat <- 470
  n <- (470+823)
  Xsample <- Xsample

  J <- 34

  X_dat <- data.frame(model.matrix(~ -1 + ., Xsample))
  num_x <- ncol(X_dat)

  # treatment ----
  xmat <- model.matrix(~ -1 + ., data.frame(X_dat=X_dat))
  alpha.X <- sqrt(0.15) # higher pretest, more likely to participate
  bx <- alpha.X / sqrt(ncol(xmat)) / ((1:ncol(xmat))^(0))

  treated_p <- n_treat / n
  ttstar <- xmat %*% bx + rnorm(n)
  tt <- cut(ttstar, quantile(ttstar, c(0, (1-treated_p), 1)), include.lowest = T)
  tt <- ifelse(as.numeric(tt)==2, 1, 0)
  table(tt)
  # Potential cluster assignment ---------------------
  wmat <- xmat

  # Observed class sizes were capped at 15 (range = 8–15, Mdn = 12).
  set.seed(12)
  bw <- sqrt(0.4) / ((1:ncol(wmat))^(0.9))

  tryk <- FALSE
  while(tryk == F) {
    K1star <- wmat %*% bw + rnorm(n)
    length( c(-6, rep(-1, 15), rep(1, 15)) )
    obs_size.k <- n_treat /J + c(-6, rep(-0.5, 16), rep(1, 6+8), rep(0,3))
    obs_prop.k <- cumsum(obs_size.k) / n_treat
    K1 <- cut(K1star, c(-Inf, quantile((K1star[tt==1]), c(obs_prop.k[1:(J-1)])), Inf), include.lowest = T)
    K1 <- as.numeric(K1)
    tryk <- min(table(K1[tt==1])) <=8 & max(table(K1[tt==1]))<=15 & min(table(K1[tt==1])) >=5
    quantile(table(K1[tt==1]))
    sum(table(K1))
  }

  ICC_Y1 <- 0.3
  sigmau1 <- sqrt( ICC_Y1 )
  sigmae1 <- sqrt( 1-ICC_Y1 )

  # cluster mean of baseline covariates
  cluster_K <- aggregate(data.frame(Xb_dat=X_dat), by = list(K1=K1), mean)
  fixed <- TRUE
  if (fixed == TRUE) {
    set.seed(12345)
    cluster_K$u1 <- (rnorm(J, 0, sigmau1))
  }

  # data so far
  id <- 1:n
  dat <- merge(data.frame(id=id, K1=K1, tt=tt, X_dat=X_dat),
               cluster_K, by = "K1")


  # Potential outcomes ----

  ratio_ve01 <- 1
  sigmae0 <- sqrt( sigmae1^2*ratio_ve01 )

  revar <- sqrt( 0.5*(sigmau1^2+sigmae1^2) + 0.5*sigmae0^2 )
  ICC_X <- ICC_W <- 0.1
  bx_control <- sqrt( revar^2*0.26/1/0.5 )
  bx_within <- sqrt( revar^2*0.15/(1-ICC_X)/0.5 )
  bx_between<- sqrt( revar^2*0.39/ICC_X/0.5 )

  e0 <- rnorm(n,0,sigmae0)
  e1 <- rnorm(n,0,sigmae1)

  gamma0 <- 0;
  Trteff.size <- 0.2

  sdpool <- sqrt(bx_control^2+bx_control^2+sigmae0^2)

  gamma1 <- Trteff.size*sdpool + gamma0

  # treatment outcome

  kmat <- cbind(
    u1 = dat$u1,
    (dat[, grep("Xb_dat", colnames(dat),value = T)] - 0)/sqrt(1),
    (dat[, grep("X_dat", colnames(dat),value = T)] - 0)/sqrt(1)
  )

  bk_treat <- 0 / ((1:ncol(kmat))^(0))
  names(bk_treat) <- colnames(kmat)

  bk_treat[grep("u1", colnames(kmat), value = T)] <- 1

  bk_treat[grep("Xb_dat", colnames(kmat), value = T)] <- (bx_between-bx_within) / c(1:length(grep("Xb_dat", colnames(kmat), value = T)))^0.9
  bk_treat[grep("X_dat", colnames(kmat), value = T)] <- bx_within / c(1:length(grep("X_dat", colnames(kmat), value = T)) )^0.9

  # homogeneous treatment effects across clusters
  if (te_k_homo == TRUE) {
    bx_control <- bx_between
    bk_treat[grep("u1", colnames(kmat), value = T)] <- 0
    e1 <- e1*0
    e0 <- e0*0
  }

  # control outcome
  cmat <- cbind(
    model.matrix(~ -1 + ., dat[, grep("X_dat", colnames(dat),value = T)])
  )
  bk_control <- 0 / ((1:ncol(cmat))^(0))
  names(bk_control) <- colnames(cmat)
  bk_control[grep("X_dat", colnames(cmat), value = T)] <- bx_control / c(1:length(grep("X_dat", colnames(cmat), value = T)))^0.9


  Y1 <- gamma1 + as.matrix(kmat) %*% bk_treat + e1
  Y0 <- gamma0 + as.matrix(cmat) %*% bk_control + e0

  # True Values
  ate.sample <- mean( Y1-Y0 )
  ate_d <- ate.sample /  sd(as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt)))

  sample.te_K <- rep(0, J)
  for(k in 1:J) { # E[Y1 - Y0 | K1=k]
    sample.te_K[k] <- mean(( Y1-Y0 )[dat$K1==k])
  }

  # Observed Data ----------------------------
  dat[1:3, ]
  dat$Y <- as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt))
  dat$K <- dat$K1 * dat$tt + 0 * (1-dat$tt)

  meanYobs <- mean(as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt)))
  sdYobs <- sd(as.numeric(Y1 * dat$tt + Y0 * (1-dat$tt)))

  datobs <- dat[, c("Y", "K", "tt",
                    grep("X_dat", colnames(dat),value = T))]

  datobs$id <- 1:nrow(datobs)


  out <- mget( ls(), envir = environment() )

  return(out)
}



OneDataSet <- function(iseed = 1, cond = 1) {
  # homogeneous ate_k --------
  gen_data_homo <- real.mimic(te_k_homo = T)

  data_homo <- gen_data_homo$datobs
  data_homo$Y <- (data_homo$Y - mean(data_homo$Y)) / sd(data_homo$Y)

  Xnames <- c(grep("X_dat", colnames(data_homo), value = T))

  suppressWarnings(
    res_homo <- cluster.specific.ate(
      data_in = data_homo,
      Xnames = Xnames,
      estimator = c("triply-robust (dml)")
    )
  )
  
  est_teK_homo <- res_homo[["triply_dml.ate_K"]]


  # heter ate_k --------
  gen_data_heter <- real.mimic(te_k_homo = F)

  data_heter <- gen_data_heter$datobs
  data_heter$Y <- (data_heter$Y - mean(data_heter$Y)) / sd(data_heter$Y)

  Xnames <- c(grep("X_dat", colnames(data_heter), value = T))

  # crossfit
  suppressWarnings(
    res_heter <- cluster.specific.ate(
      data_in = data_heter,
      Xnames = Xnames,
      estimator = c("triply-robust (dml)")
    )
    
  )
  
  est_teK_heter <- res_heter[["triply_dml.ate_K"]]

  library(tidyverse)
  # Plot distribution ----------
  fig.theme <-
    theme_bw() +
    theme(panel.grid.minor = element_line(linewidth = 0),
          panel.grid.major.x = element_line(linewidth = 0),
          panel.grid.major.y = element_line(linewidth = 0.5, lineend = "round", color = "grey", linetype = "longdash"),
          # axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.spacing.x = unit(0.1, "mm"),
          legend.key.height = unit(15, "mm"),
          legend.key.size = unit(10, "mm"))


  # Histogram
  histdata <- rbind(data.frame(est_teK_homo, cluster_specific_effect = "truly homogeneous"),
                    data.frame(est_teK_heter, cluster_specific_effect = "truly heterogenous"))
  sumdata <- histdata %>% group_by(cluster_specific_effect) %>% summarise(
    SD_ate_k = sd(ate_k),
    min = min(ate_k), max = max(ate_k), median = median(ate_k),
    negeff = mean(ate_k < 0), sum(ate_k < 0)
  )
  hist_ate.k <- ggplot(histdata, aes(x = ate_k)) +
    geom_histogram(aes(y = after_stat(ncount)), bins = 20, fill="grey") +
    geom_density(aes(y = after_stat(scaled)), linewidth = 0.5 ) +
    scale_x_continuous("Effect estimates across teacheres/classes (clusters) in the treatment arm (ATE_k)" #, limits = c(-2,3)
    ) +
    scale_y_continuous("Proportion") +
    geom_text(data = sumdata,
              mapping = aes(
                label = paste0("SD of ATE_k = ", round(SD_ate_k,2)),
                x = -0.9, y = 1.1
              ), size = 4.5) +
    facet_grid(. ~ cluster_specific_effect, scales = "fixed", labeller = "label_both") +
    fig.theme
  hist_ate.k

  # Caterpillar plot
  pdat <- rbind(
    data.frame(est_teK_homo[order(est_teK_homo$ate_k), ], cluster_specific_effect = "truly homogeneous",
               cluster = 1:nrow(est_teK_homo)),
    data.frame(est_teK_heter[order(est_teK_heter$ate_k), ], cluster_specific_effect = "truly heterogenous",
               cluster = 1:nrow(est_teK_heter))
  )

  ciplot_ate.k <- ggplot( data = pdat ) +
    geom_point(aes(x = cluster, y = ate_k), shape = 0, size=3) +
    geom_segment(aes(x = cluster, xend = cluster, y = ate_k_ci1, yend = ate_k_ci2), linewidth = 1, inherit.aes = F) +
    scale_x_discrete("Teacher/Class (i.e., cluster) in the treatment arm") +
    scale_y_continuous("Cluster-specific effects: Estimates and 95% CI") +
    facet_grid(. ~ cluster_specific_effect, scales = "fixed", labeller = "label_both") +
    fig.theme

  ciplot_ate.k

  example_out <- list(
    dataset_homo = data_homo,
    dataset_heter = data_heter,
    est_teK_homo = est_teK_homo,
    est_teK_heter = est_teK_heter,
    hist_ate.k = hist_ate.k,
    ciplot_ate.k = ciplot_ate.k
  )

  return(example_out)
}

save(example_out, file = "~/Library/CloudStorage/Box-Box/Labs/Causal/PND/R_github/example_out.RData")
