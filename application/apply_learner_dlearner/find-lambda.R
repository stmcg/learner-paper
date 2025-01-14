################################################################################
## Loading packages and libraries
################################################################################
library('MASS')
library('learner')

load('../dat.RData')

################################################################################
## Global parameters
################################################################################

set.seed(1234)

n_lambdas <- 5
lambda_1_vals <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -2.5, to = 2.5, length.out = n_lambdas)
n_cores <- 13
max_iter <- 100
c <- 0.025

################################################################################
## Applying LEARNER
################################################################################

fit <- cv.learner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen, 
                  lambda_1_all = lambda_1_vals, 
                  lambda_2_all = lambda_2_vals, 
                  step_size = c, n_cores = n_cores, 
                  control = list(max_iter = max_iter))

out_small_lambda <- list(res = fit$mse_all, 
                         lambda1_min = fit$lambda_1_min, 
                         lambda2_min = fit$lambda_2_min)

save(out_small_lambda, file = 'lambda_small.RData')
