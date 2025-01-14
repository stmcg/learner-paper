################################################################################
## Loading packages and libraries
################################################################################
rm(list = ls())

library('MASS')
library('learner')

################################################################################
## Global parameters
################################################################################

load('../dat.RData')
max_iter <- 100
c <- 0.025

load('lambda_small.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min

fit <- try(learner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen, 
                   lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
                   step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)

save(fit, file = 'learner.RData')
