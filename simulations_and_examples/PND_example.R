# loading packages ----------------------
library(tidyverse)
library(labelled)
library(purrr)
library(origami)
library(mvtnorm)
library(glue)
library(SuperLearner)

# loading the method package 
# devtools::install_github(repo = "xliu12/PND.heter", subdir = "PND.heter.cluster")
# library(PND.heter.cluster)
# or download the package and then
devtools::load_all("../PND.heter.cluster")

# loading the data
load("../simulations_and_examples/data_example.RData")
# the variables
ttname <- grep("^tt", colnames(data_in), value = T)
Kname <- grep("^K", colnames(data_in), value = T)
Yname <- grep("^Y_TICINSUPMN", colnames(data_in), value = T)
Xnames <- grep("^X_", colnames(data_in), value = T)

Yfamily <- ifelse(length(unique(data_in[[Yname]])) > 2, "gaussian", "binomial")
## run -------------
set.seed(12345)
crossfit_res <- atekCl(data_in = data_in,
                       ttname = ttname, Kname = Kname, Yname = Yname,
                       Xnames = Xnames,
                       Yfamily = "gaussian",
                       learners_tt = c("SL.nnet", "SL.ranger"),
                       learners_k = c("SL.nnet.modified"),
                       learners_y = c("SL.nnet", "SL.ranger"),
                       sensitivity = "small_to_medium")


## covariate balance -------
cov_balance <- balance(data_in = data_in, 
                       atekCl_results = crossfit_res, 
                       Yfake_name = Xnames, 
                       ttname = ttname, Kname = Kname) 

summ_cov_balance <- cov_balance %>% 
  group_by(covariate) %>% 
  summarise(across(contains("smd"), list(
    max_abs = ~max(abs(.)), median_abs = ~quantile(abs(.), 0.5)
  )))

covariate_name <- c("X_agreeableness", "X_conscientiousness", "X_extraversion", "X_class_poverty", "X_teaching_experience", "X_headstart", "X_public", "X_emotional_support","X_instructional_support", "X_organizational_support", "X_age", "X_gender", "X_income/needs_ratio", "X_self_efficacy", "X_race_Black", "X_race_Hispanic","X_race_White", "X_yrs_education", "X_parent_edu","X_socioeconomic_status")
summ_cov_balance$covariate_name <- covariate_name



cov_balance1 <- cov_balance %>% 
  full_join(summ_cov_balance[, c("covariate", "covariate_name")], by = "covariate") 
## plot x ------
cov_balance_long <- cov_balance1 %>% 
  pivot_longer(cols = c(contains("smd")), names_to = "if_adj", values_to = "smd") %>% 
  mutate(if_adj = factor(if_adj, levels=c("smd", "unadj_smd"), labels=c("After adjustment", "Unadjusted"))) %>% 
  select(covariate,covariate_name, everything())

bal.plot <- cov_balance_long %>% 
  ggplot(aes(x = covariate_name, y = abs(smd))) +
  # geom_point() +
  geom_boxplot(aes(fill = if_adj), position = "dodge", outliers = F) +
  geom_hline(yintercept = 0.1) +
  labs(x="Covariates", y="|Mean Diff| / SD") +
  scale_fill_discrete("Covariate balance between each cluster of the treatment arm vs. the control arm ") +
  scale_y_continuous(breaks = c(0, seq(0.1, 2 ,0.3))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, size = 8, hjust = 0.65),
    axis.text = element_text(family = "sans"),
    axis.title = element_text(size = 15, family = "sans"),
    legend.position = "top", 
    legend.text = element_text(size = 15, family = "sans"),
    legend.title = element_text(size = 15, family = "sans"),
    legend.direction = "vertical",
    legend.box = "vertical",
    panel.grid.minor.y = element_blank()
  )
bal.plot


## cluster-specific treatment results ----

ate_K_eg <- data.frame(crossfit_res$ate_K, 
                        cluster_K = unique(data_in[[Kname]][data_in[[ttname]]==1]) ) 

sens_ate_K_eg <- purrr::map(1:length(crossfit_res$sens_results), 
                             function(l=1) {
                               resl1 <- data.frame(sens_smd=as.numeric(names(crossfit_res$sens_results)[l]), 
                                                   cluster_K = unique(data_in[[Kname]][data_in[[ttname]]==1]),
                                                   crossfit_res$sens_results[[l]])
                               resl1
                             } ) %>% 
  reduce(bind_rows)


## main plot ----
ciplot <- ate_K_eg %>% 
  arrange(ate_k) %>% 
  mutate(cluster = factor(order(ate_k)), 
         ci1=boot_ci1, ci2=boot_ci2) %>% 
  ggplot(aes(y = ate_k, x= cluster)) +
  geom_pointrange(aes(ymin = ci1, ymax = ci2)) + 
  geom_errorbar(aes(ymin = ci1, ymax = ci2)) +
  # geom_density(aes(y = after_stat(scaled)), linewidth = 0.5 ) +
  scale_x_discrete("Coaches (i.e., clusters) in the treatment arm") + 
  scale_y_continuous("Estimated treatment effect") +
  geom_hline( yintercept = 0, color="brown", linetype="solid", linewidth=1) +
  # geom_text(data = sumdata,
  #           mapping = aes(
  #             label = paste0("SD of ATE_k = ", round(SD_ate_k,2)), 
  #             x = -0.9, y = 1.1
  #           ), size = 4.5) +
  # facet_grid(. ~ cluster_specific_effect, scales = "fixed", labeller = "label_both") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
ciplot



## sensitivity plot ---------
J <- length(unique(ate_K_eg$cluster_K))

sensPlot_eg <- sens_ate_K_eg %>% 
  # mutate(cluster = cluster_K) %>% 
  mutate(
    cluster = factor(cluster_K, levels=order(ate_K_eg$ate_k), labels=glue("ATE_cluster_{1:J}")),
    ci1 = boot_ci1, ci2 = boot_ci2 ) %>% 
  ggplot(aes(x = sens_smd, y = SenSate_k, ymin=ci1, ymax=ci2)) +
  facet_wrap(~ cluster, scales = "fixed") +
  geom_pointrange(aes(ymin = ci1, ymax = ci2)) + 
  geom_errorbar(aes(ymin = ci1, ymax = ci2)) +
  # geom_density(aes(y = after_stat(scaled)), linewidth = 0.5 ) +
  scale_x_continuous("Sensitivity parameter: Mean difference standardized"
                     , breaks = unique(sens_ate_K_eg2$sens_smd)) + 
  scale_y_continuous("Estimated treatment effect") +
  geom_hline( yintercept = 0, color="brown", linetype="solid", linewidth=1) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )
sensPlot_eg
