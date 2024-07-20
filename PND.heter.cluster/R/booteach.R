

# individual bootstrap scores


fun.each.boot <- function(cv_components, boot_inds, data_in) {
  
  boot_components <- cv_components[boot_inds, ]
  
  boot_ate_K <- lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
    # bootstrap scores for estimating ate_k
    ymean.kt1 <- boot_components[, paste0("ymean.kt1_dr_valid.V", k)]
    ymean.t0 <- boot_components[, paste0("ymean.t0_dr_valid.V", k)]
    pk <- boot_components[, paste0("pk_dr_valid.V", k)]
    
    boot_ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)
    
    return(c(boot_ate_k = boot_ate_k))
  })
  boot_ate_K <- data.frame(do.call(rbind, boot_ate_K))
  each_boot <- c(boot_ate_k = t(boot_ate_K$boot_ate_k))
  
  # # if do linear combination
  # if (length(combination)>0) {
  #   if (combination[1]=="cluster_average") {
  #     num_clus <- length( unique(data_in$K[data_in$tt==1]) )
  #     comb <- rep(1/num_clus, num_clus)
  #   }
  #   if (combination[1]!="cluster_average") {
  #     comb <- combination
  #   }
  #   boot_ate_comb <- sum(comb * boot_ate_K$boot_ate_k)
  #   each_boot <- c(boot_ate_k = t(boot_ate_K$boot_ate_k), boot_ate_comb=boot_ate_comb)
  #   
  # }
  
  
  
  return(each_boot)
}

fun.each.bootsens <-  function(cv_components, boot_inds, dk, data_in) {
  
  boot_components <- cv_components[boot_inds, ]
  
  boot_sens_ate_K <- do.call(rbind, lapply(unique(data_in$K[data_in$tt==1]), FUN = function(k=1) {
    
    ymean.kt1 <- boot_components[, paste0("ymean.kt1_dr_valid.V", k)]
    ymean.t0 <- boot_components[, paste0("ymean.t0_dr_valid.V", k)]
    pk <- boot_components[, paste0("pk_dr_valid.V", k)]
    pkpk <- boot_components[, paste0("pkpk_dr_valid.V", k)]
    
    # point estimate
    ate_k <- mean(ymean.kt1 - ymean.t0) / mean(pk)
    psi_sens_k <- mean(pkpk)/mean(pk)
    SenSate_k <- ate_k - dk + dk*psi_sens_k
    
    
    c(SenSate_k = SenSate_k)
  }))
  
  each_bootsens <- c(boot_sens = as.numeric(boot_sens_ate_K))
  
  each_bootsens
}


