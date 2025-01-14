## Loading data
rm(list = ls())
load('../dat.RData')

library('softImpute')

## Estimating singular spaces in source population
svd1 <- svd(theta_hat_1, nu = r_chosen, nv = r_chosen)
PU_1 <- svd1$u %*% t(svd1$u)
PV_1 <- svd1$v %*% t(svd1$v)

## Set 1
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[1]]] <- NA
fit <- softImpute(x = theta_hat_0_red, lambda = 0, rank.max = r_chosen, maxit = 200)
theta_hat_0_imputed <- complete(theta_hat_0_red, fit)
hard_set1 <- theta_hat_0_imputed[heldout_set[[1]]]
dlearner_set1 <- (PU_1 %*% theta_hat_0_imputed %*% PV_1)[heldout_set[[1]]]

save(hard_set1, file = 'hard_set1.RData')
save(dlearner_set1, file = 'dlearner_set1.RData')

## Set 2
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[2]]] <- NA
fit <- softImpute(x = theta_hat_0_red, lambda = 0, rank.max = r_chosen, maxit = 200)
theta_hat_0_imputed <- complete(theta_hat_0_red, fit)
hard_set2 <- theta_hat_0_imputed[heldout_set[[2]]]
dlearner_set2 <- (PU_1 %*% theta_hat_0_imputed %*% PV_1)[heldout_set[[2]]]

save(hard_set2, file = 'hard_set2.RData')
save(dlearner_set2, file = 'dlearner_set2.RData')

## Set 3
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[3]]] <- NA
fit <- softImpute(x = theta_hat_0_red, lambda = 0, rank.max = r_chosen, maxit = 200)
theta_hat_0_imputed <- complete(theta_hat_0_red, fit)
hard_set3 <- theta_hat_0_imputed[heldout_set[[3]]]
dlearner_set3 <- (PU_1 %*% theta_hat_0_imputed %*% PV_1)[heldout_set[[3]]]

save(hard_set3, file = 'hard_set3.RData')
save(dlearner_set3, file = 'dlearner_set3.RData')

## Set 4
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[4]]] <- NA
fit <- softImpute(x = theta_hat_0_red, lambda = 0, rank.max = r_chosen, maxit = 200)
theta_hat_0_imputed <- complete(theta_hat_0_red, fit)
hard_set4 <- theta_hat_0_imputed[heldout_set[[4]]]
dlearner_set4 <- (PU_1 %*% theta_hat_0_imputed %*% PV_1)[heldout_set[[4]]]

save(hard_set4, file = 'hard_set4.RData')
save(dlearner_set4, file = 'dlearner_set4.RData')

## Set 5
theta_hat_0_red <- theta_hat_0
theta_hat_0_red[heldout_set[[5]]] <- NA
fit <- softImpute(x = theta_hat_0_red, lambda = 0, rank.max = r_chosen, maxit = 200)
theta_hat_0_imputed <- complete(theta_hat_0_red, fit)
hard_set5 <- theta_hat_0_imputed[heldout_set[[5]]]
dlearner_set5 <- (PU_1 %*% theta_hat_0_imputed %*% PV_1)[heldout_set[[5]]]

save(hard_set5, file = 'hard_set5.RData')
save(dlearner_set5, file = 'dlearner_set5.RData')
