rm(list = ls())
source('helper.R')
source('helper-unequal-rank.R')
n_iter <- 50

# Simulation scenarios with unequal ranks
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c_r0_4 <- 0.0075; c_r0_8 <- 0.0225
max_iter <- 75

## Setting with r_0 = 4 and r_1 = 8
set.seed(12345)
title_help <- 'unequal-r0=4-r1=8'
res_4_8 <- run_all_scenarios_unequal_rank(var_0_all = var_0_all, var_1_all = var_1_all, 
                                          r_0 = 4, r_1 = 8, r_all = c(4, 8), independent_noise = TRUE)

## Setting with r_0 = 8 and r_1 = 4
set.seed(12345)
title_help <- 'unequal-r0=8-r1=4'
res_8_4 <- run_all_scenarios_unequal_rank(var_0_all = var_0_all, var_1_all = var_1_all, 
                                          r_0 = 8, r_1 = 4, r_all = c(4, 8), independent_noise = TRUE)

save.image(file = 'simres-unequal-rank.RData')
