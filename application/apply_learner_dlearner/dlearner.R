load('../dat.RData')

svd1 <- svd(theta_hat_1, nu = r_chosen, nv = r_chosen)

PU_1 <- svd1$u %*% t(svd1$u)
PV_1 <- svd1$v %*% t(svd1$v)
theta_hat_0_dlearner <- PU_1 %*% theta_hat_0 %*% PV_1

colnames(theta_hat_0_dlearner) <- colnames(theta_hat_1)
rownames(theta_hat_0_dlearner) <- rownames(theta_hat_1)

save(theta_hat_0_dlearner, file = 'DLEARNER-estimate-BBJ.RData')