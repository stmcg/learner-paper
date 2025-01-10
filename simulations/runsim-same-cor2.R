rm(list = ls())
source('helper.R')

# Simulation scenarios with same latent spaces
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c <- 0.035
max_iter <- 75
title_help <- 'same-rho=0.25'

set.seed(12345)
rho <- 0.25
res_same_cor2 <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 4, 
                                   same_latent_space = TRUE, independent_noise = FALSE)

save.image(file = 'simres-same-cor2.RData')
