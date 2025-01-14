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

## Set 1
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[1]]] <- NA
load('lambda_small_set1.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min
print(lambda1_min); print(lambda2_min)
fit <- try(
  learner(Y_source = theta_hat_1, Y_target = theta_hat_0_red, r = r_chosen, 
          lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
          step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)
learner_set1 <- fit$learner_estimate[heldout_set[[1]]]
save(learner_set1, file = 'learner_set1.RData')

## Set 2
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[2]]] <- NA
load('lambda_small_set2.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min
print(lambda1_min); print(lambda2_min)
fit <- try(
  learner(Y_source = theta_hat_1, Y_target = theta_hat_0_red, r = r_chosen, 
          lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
          step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)
learner_set2 <- fit$learner_estimate[heldout_set[[2]]]
save(learner_set2, file = 'learner_set2.RData')

## Set 3
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[3]]] <- NA
load('lambda_small_set3.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min
print(lambda1_min); print(lambda2_min)
fit <- try(
  learner(Y_source = theta_hat_1, Y_target = theta_hat_0_red, r = r_chosen, 
          lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
          step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)
learner_set3 <- fit$learner_estimate[heldout_set[[3]]]
save(learner_set3, file = 'learner_set3.RData')

## Set 4
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[4]]] <- NA
load('lambda_small_set4.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min
print(lambda1_min); print(lambda2_min)
fit <- try(
  learner(Y_source = theta_hat_1, Y_target = theta_hat_0_red, r = r_chosen, 
          lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
          step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)
learner_set4 <- fit$learner_estimate[heldout_set[[4]]]
save(learner_set4, file = 'learner_set4.RData')

## Set 5
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[5]]] <- NA
load('lambda_small_set5.RData')
lambda1_min <- out_small_lambda$lambda1_min
lambda2_min <- out_small_lambda$lambda2_min
print(lambda1_min); print(lambda2_min)
fit <- try(
  learner(Y_source = theta_hat_1, Y_target = theta_hat_0_red, r = r_chosen, 
          lambda_1 = lambda1_min, lambda_2 = lambda2_min, 
          step_size = c, control = list(max_iter = max_iter)))
print(fit$objective_values)
learner_set5 <- fit$learner_estimate[heldout_set[[5]]]
save(learner_set5, file = 'learner_set5.RData')
