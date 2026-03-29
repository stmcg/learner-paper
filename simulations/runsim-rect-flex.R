rm(list = ls())
source('helper.R')
source('helper-rect-flex.R')

n_lambdas <- 5
max_iter <- 75

# ------------------------------------------------------------------------------
# High similarity: same latent space (identical U and V across sources)
# Lambda / step analogous to runsim-same.R
# ------------------------------------------------------------------------------
lambda_1_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c <- 0.0035
title_help <- 'rect-flex-high'

set.seed(1234)
res_rect_high <- run_all_scenarios_rect_flex(var_0_all = var_0_all, var_1_all = var_1_all, r = 4,
                                             same_latent_space = TRUE, dif_factor = 4)
save(res_rect_high, 
     file = 'simres-rect-flex-high.RData', compress = 'xz')

# ------------------------------------------------------------------------------
# Moderate similarity: different latent spaces, dif_factor = 4 (both U and V perturbed)
# Lambda / step analogous to runsim-dif-moderate.R
# ------------------------------------------------------------------------------
lambda_1_vals <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -2, to = 1, length.out = n_lambdas)
c <- 0.035
title_help <- 'rect-flex-moderate'

set.seed(1234)
res_rect_moderate <- run_all_scenarios_rect_flex(var_0_all = var_0_all, var_1_all = var_1_all, r = 4,
                                                 same_latent_space = FALSE, dif_factor = 4)

save(res_rect_moderate, 
     file = 'simres-rect-flex-moderate.RData', compress = 'xz')
