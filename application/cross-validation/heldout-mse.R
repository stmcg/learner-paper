## Loading data
rm(list = ls())
load('../dat.RData')

## Loading LEARNER predictions
load('learner_set1.RData')
load('learner_set2.RData')
load('learner_set3.RData')
load('learner_set4.RData')
load('learner_set5.RData')

## Loading D-LEARNER predictions
load('dlearner_set1.RData')
load('dlearner_set2.RData')
load('dlearner_set3.RData')
load('dlearner_set4.RData')
load('dlearner_set5.RData')

## Loading Target-Only SVD predictions
load('hard_set1.RData')
load('hard_set2.RData')
load('hard_set3.RData')
load('hard_set4.RData')
load('hard_set5.RData')

## Getting Source-Only SVD estimates
svd1 <- svd(theta_hat_1, nu = r_chosen, nv = r_chosen)
theta_hat_1_svd <- svd1$u %*% diag(svd1$d[1:r_chosen], nrow = r_chosen, ncol = r_chosen) %*% t(svd1$v)
svd1_set1 <- theta_hat_1_svd[heldout_set[[1]]]
svd1_set2 <- theta_hat_1_svd[heldout_set[[2]]]
svd1_set3 <- theta_hat_1_svd[heldout_set[[3]]]
svd1_set4 <- theta_hat_1_svd[heldout_set[[4]]]
svd1_set5 <- theta_hat_1_svd[heldout_set[[5]]]

## Getting true values
true_set1 <- theta_hat_0[heldout_set[[1]]]
true_set2 <- theta_hat_0[heldout_set[[2]]]
true_set3 <- theta_hat_0[heldout_set[[3]]]
true_set4 <- theta_hat_0[heldout_set[[4]]]
true_set5 <- theta_hat_0[heldout_set[[5]]]

## Getting performance
round(mean((svd1_set1 - true_set1)^2), 2)
round(mean((hard_set1 - true_set1)^2), 2)
round(mean((dlearner_set1 - true_set1)^2), 2)
round(mean((learner_set1 - true_set1)^2), 2)

round(mean((svd1_set2 - true_set2)^2), 2)
round(mean((hard_set2 - true_set2)^2), 2)
round(mean((dlearner_set2 - true_set2)^2), 2)
round(mean((learner_set2 - true_set2)^2), 2)

round(mean((svd1_set3 - true_set3)^2), 2)
round(mean((hard_set3 - true_set3)^2), 2)
round(mean((dlearner_set3 - true_set3)^2), 2)
round(mean((learner_set3 - true_set3)^2), 2)

round(mean((svd1_set4 - true_set4)^2), 2)
round(mean((hard_set4 - true_set4)^2), 2)
round(mean((dlearner_set4 - true_set4)^2), 2)
round(mean((learner_set4 - true_set4)^2), 2)

round(mean((svd1_set5 - true_set5)^2), 2)
round(mean((hard_set5 - true_set5)^2), 2)
round(mean((dlearner_set5 - true_set5)^2), 2)
round(mean((learner_set5 - true_set5)^2), 2)
