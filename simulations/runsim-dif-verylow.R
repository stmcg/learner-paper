rm(list = ls())
source('helper.R')

# Simulation scenarios with different latent spaces
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -2, to = 1, length.out = n_lambdas)
c <- 0.07
max_iter <- 75
title_help <- 'dif-verylow'

set.seed(1234)
res_1_dif_verylow <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 4, 
                                       same_latent_space = FALSE, dif_factor = 2, independent_noise = TRUE)
res_2_dif_verylow <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 8, 
                                       same_latent_space = FALSE, dif_factor = 2, independent_noise = TRUE)

save.image(file = 'simres-dif-verylow.RData')