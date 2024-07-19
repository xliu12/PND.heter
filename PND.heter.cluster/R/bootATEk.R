


# bootstrap-based with cv_components
boot.atek.cl <- function(seedb=123,
    data_in,
    cv_components,
    Xnames, ttname, Kname, Yname,
    Yfamily
) {
  # cv_components <- crossfit_res$cv_components
  # colnames(cv_components)
  set.seed(seedb)
  boot_inds <- sample(1:nrow(data_in), nrow(data_in),replace = T)
  X_boot <- data_in[boot_inds, Xnames]
  # (X,tt,K1) independent of each other
  n <- nrow(data_in)
  J <- length(unique(data_in$K[data_in$tt==1]))
  pt_boot <- cv_components$t1.x_pred_valid[boot_inds]
  # tt_boot <- rbinom(n, 1, pt_boot)
  tt_boot <- data_in[boot_inds, ttname, drop=TRUE]

  pk.multi_boot <- cv_components[boot_inds, glue("pk.multi_pred_valid.V{1:J}")]
  # mChoices_boot <- t(apply(pk.multi_boot, 1, rmultinom, n=1, size=1))
  # K1_boot <- apply(mChoices_boot, 1, function(x) {which(x==1)})
  # K_boot <- K1_boot*1*(tt_boot==1) + 0*(tt_boot==0)
  K_boot <- data_in[boot_inds, Kname, drop=TRUE]

  nkboot_min <- min(table(K_boot[tt_boot==1]))
  while(nkboot_min < 2) {
    tt_boot <- rbinom(n, 1, pt_boot)
    nkboot_min <- min(table(K_boot[tt_boot==1]))
    # print("aa")
  }


  df_boot <- data.frame(boot_inds=boot_inds, id_boot=1:n, tt_boot=tt_boot, K_boot=K_boot)
  df_boot$Y_boot <- NA
  # Y|K,X,tt=1
  k <-1
  cv_components$y1.Kobs_pred <- NA
  for(k in 1:J) {
    cv_components$y1.Kobs_pred[data_in[[Kname]]==k] <-
      cv_components[data_in[[Kname]]==k, glue("y.kt1_pred_valid.V{k}")]
  }
  obs_resi1 <- data_in[[Yname]][data_in[[ttname]]==1] - cv_components[data_in[[ttname]]==1, glue("y1.Kobs_pred")]

  # reorder cv_components
  cv_components_rep <- cv_components #[df_boot$boot_inds, ]

  for(k in 1:J) {
    if(Yfamily=="gaussian") {
      # obs_resi1k <- data_in[[Yname]][data_in[[Kname]]==k] - cv_components[data_in[[Kname]]==k, glue("y.kt1_pred_valid.V{k}")]
      obs_resi1k <- obs_resi1

      y.kt1_boot <- sample(obs_resi1k, sum(K_boot==k), replace = TRUE) +
        # cv_components_rep[K_boot==k, glue("y.kt1_pred_valid.V{k}")]
        sample(cv_components_rep[, glue("y.kt1_pred_valid.V{k}")], sum(K_boot==k), replace = TRUE)

      # y.kt1_boot <- rnorm(sum(K_boot==k), 0, sd = sd(obs_resi1k)) +
      #   cv_components[K_boot==k, glue("y.kt1_pred_valid.V{k}")]
    }
    if(Yfamily=="binomial") {
      y.kt1_boot <- rbinom(sum(K_boot==k), 1, cv_components_rep[K_boot==k, glue("y.kt1_pred_valid.V{k}")])
    }

    df_boot$Y_boot[df_boot$K_boot==k] <- y.kt1_boot
  }
  # Y|X,tt=0
  if(Yfamily=="gaussian") {
    obs_resi0 <- data_in[[Yname]][data_in[[ttname]]==0] - cv_components[data_in[[ttname]]==0, glue("y.t0xw_pred_valid")]

    y.t0_boot <- sample(obs_resi0, sum(tt_boot==0), replace = TRUE) +
      # cv_components_rep[tt_boot==0, glue("y.t0xw_pred_valid")]
      sample(cv_components_rep[, glue("y.t0xw_pred_valid")], sum(tt_boot==0), replace = TRUE)

    # y.t0_boot <- rnorm(sum(tt_boot==0), 0, sd = sd(obs_resi0)) +
    #   cv_components[tt_boot==0, glue("y.t0xw_pred_valid")]
  }
  if(Yfamily=="binomial") {
    y.t0_boot <- rbinom(sum(tt_boot==0), 1, cv_components_rep[tt_boot==0, glue("y.t0xw_pred_valid")])
  }
  df_boot$Y_boot[df_boot$tt_boot==0] <- y.t0_boot

  # estimator with boot
  wt1_boot <- 1*(tt_boot==1)/(pt_boot)
  wt0_boot <- 1*(tt_boot==0)/(1-pt_boot)

  pk.dr_boot <- map(1:J, function(k=1) {
    pk_dr_boot <- wt1_boot*1*(tt_boot==1)*(1*(K_boot==k) - pk.multi_boot[[k]]) + pk.multi_boot[, k]

    pk_dr_boot
  })

  y.kt1.dr_boot <- map(1:J, function(k=1) {
    y.kt1_pred_boot <- cv_components[df_boot$boot_inds, glue("y.kt1_pred_valid.V{k}")]

    y.kt1_dr_boot <- wt1_boot*1*(tt_boot==1)*(K_boot==k)*(df_boot$Y_boot - y.kt1_pred_boot) +
      #y.kt1_pred_boot * cv_components[df_boot$boot_inds, glue("pk_dr_valid.V{k}")]
      y.kt1_pred_boot * pk.dr_boot[[k]]

    y.kt1_dr_boot
  })

  y.kt0.dr_boot <- map(1:J, function(k=1) {
    # y.t0_pred_boot <- cv_components[df_boot$boot_inds, glue("y.t0xw_pred_valid")]
    # y.t0_dr_boot <- wt0_boot*1*(tt_boot==0)*pk.multi_boot[[k]]*(df_boot$Y_boot - y.t0_pred_boot) + y.t0_pred_boot * pk.dr_boot[[k]]

    y.t0_dr_boot <- cv_components[df_boot$boot_inds, glue("ymean.t0_dr_valid.V{k}")]
    y.t0_dr_boot
  })

  ate.k_boot <- map_dbl(1:J, function(k=1) {
    pkdr_boot <- pk.dr_boot[[k]]
    # pkdr_boot <- cv_components[df_boot$boot_inds, glue("pk_dr_valid.V{k}")]
    atek_boot <- mean(y.kt1.dr_boot[[k]] - y.kt0.dr_boot[[k]]) / mean(pkdr_boot)
    atek_boot
  })

  ate.k_boot
}

atekCl <- function(nboot =0,
                   cv_folds = 4L,
                   data_in,
                   ttname = "tt", Kname = "K", Yname = "Y",
                   Xnames,
                   Fit = "mlr", # Fit = "mlr"
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

  set.seed(12345)
  crossfit_res <- cluster.specific.ate(
    cv_folds = cv_folds,
    data_in = data_in,
    ttname = ttname, Kname = Kname, Yname = Yname,
    Xnames = Xnames,
    Fit = Fit, # Fit = "mlr"
    omit.tt = omit.tt,
    omit.k = omit.k,
    y1model_lme = y1model_lme,
    Yfamily = Yfamily,
    learners_tt = learners_tt,
    learners_k = learners_k,
    learners_y = learners_y,
    combination = combination,
    sensitivity = sensitivity
  )

  if (nboot>0) {
    cv_components <- crossfit_res$cv_components

    set.seed(12333)
    seeds_boot <- sample(1:1e4, nboot, replace = F)
    boots <- sapply(seeds_boot, boot.atek.cl,
                    data_in = data_in,
                    cv_components = cv_components,
                    Xnames = Xnames,
                    ttname = ttname,
                    Kname = Kname,
                    Yname = Yname,
                    Yfamily = Yfamily)
    crossfit_res$ate_K$ate_k_var_bootcl <- apply(boots, 1, var, na.rm=TRUE)
    crossfit_res$ate_K$ate_k_ci1_bootcl <- apply(boots, 1, quantile, 0.025, na.rm=TRUE)
    crossfit_res$ate_K$ate_k_ci2_bootcl <- apply(boots, 1, quantile, 0.975, na.rm=TRUE)
  }

  return(crossfit_res)

}

