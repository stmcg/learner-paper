rm(list = ls())
source('helper.R')
source('helper-ext.R')

# Simulation scenarios with same latent spaces
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c <- 0.0075
max_iter <- 200
title_help <- 'same-rho=0.5-ext'

set.seed(12345)
rho <- 0.5
res_same_cor3_ext <- run_all_scenarios(var_0_all = var_0_all, var_1_all = var_1_all, r = 4, 
                                       same_latent_space = TRUE, independent_noise = FALSE)

save.image(file = 'simres-same-cor3-ext.RData')
