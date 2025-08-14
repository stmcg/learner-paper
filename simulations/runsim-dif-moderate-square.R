rm(list = ls())
source('helper.R')

p <- 500; q <- 500
max_rank <- 16

# Simulation scenarios with different latent spaces
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = 1, to = 3, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -6, to = 0, length.out = n_lambdas)
c <- 0.035
max_iter <- 75
title_help <- 'dif-moderate-square'

set.seed(1234)
res_1_dif_moderate_square <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 4, 
                                        same_latent_space = FALSE, independent_noise = TRUE)
res_2_dif_moderate_square <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 8,
                                        same_latent_space = FALSE, independent_noise = TRUE)

save.image(file = 'simres-dif-moderate-square.RData')