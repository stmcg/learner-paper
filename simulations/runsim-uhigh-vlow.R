rm(list = ls())
source('helper.R')
source('helper-uhigh-vlow.R')

p <- 500; q <- 500
max_rank <- 16

# Simulation scenarios: U high similarity (shared U), V similarity controlled by dif_factor_v
# (see set_theta_Ushared_Vdiff in helper-uhigh-vlow.R: threshold for V perturbation is 1/(dif_factor_v*sqrt(q))).
# Larger dif_factor_v => smaller perturbation => higher V similarity.

r <- 4
n_lambdas <- 5
max_iter <- 75

# ------------------------------------------------------------------------------
# Scenario 1: U high / V low  (dif_factor_v = 2; strongest V perturbation in this family)
# ------------------------------------------------------------------------------
dif_factor_v <- 2
lambda_1_vals <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -3, to = 1, length.out = n_lambdas)
c <- 0.07
title_help <- 'uhigh-vlow'

set.seed(1234)
res_uhigh_vlow <- run_all_scenarios_uhigh_vlow_square(
  var_0_all = var_0_all, var_1_all = var_1_all,
  r = r, dif_factor_v = dif_factor_v
)

save(res_uhigh_vlow, file = 'simres-uhigh-vlow.RData', compress = 'xz')

# ------------------------------------------------------------------------------
# Scenario 2: U high / V moderate  (dif_factor_v = 4; half the V perturbation vs vlow)
# ------------------------------------------------------------------------------
dif_factor_v <- 4
lambda_1_vals <- 10^seq(from = 1, to = 3, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -6, to = 0, length.out = n_lambdas)
c <- 0.035
title_help <- 'uhigh-vmoderate'

set.seed(1234)
res_uhigh_vmoderate <- run_all_scenarios_uhigh_vlow_square(
  var_0_all = var_0_all, var_1_all = var_1_all,
  r = r, dif_factor_v = dif_factor_v
)

save(res_uhigh_vmoderate, file = 'simres-uhigh-vmoderate.RData', compress = 'xz')
