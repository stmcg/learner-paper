# Performs the simulation study for one data generation scenario, where LEARNER uses an external data set
run_1_scenario <- function(var_0, var_1, r, thetas, independent_noise){
  res <- foreach(i = 1:n_iter) %dorng% {
    ## Data Generation
    theta_hat_0 <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
    theta_hat_1 <- get_theta_hat(thetas$theta_1, var_val = var_1, independent_noise = independent_noise)
    theta_hat_0_external <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
    
    ## Truncated SVD
    svd0 <- svd(theta_hat_0, nu = r, nv = r)
    tsvd_est <- svd0$u %*% diag(svd0$d[1:r], nrow = r, ncol = r) %*% t(svd0$v)
    res_tsvd = norm(tsvd_est - thetas$theta_0, type = 'F')
    
    ## D-LEARNER
    r_chosen <- max(adaptiveHardThresholding(Y = theta_hat_1, k = max_rank)$r, 1)
    dlearner_est <- dlearner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen)$dlearner_estimate
    res_dlearner = norm(dlearner_est - thetas$theta_0, type = 'F')
    
    ## LEARNER
    norms_training <- matrix(NA, nrow = n_lambdas, ncol = n_lambdas)
    for (l1_ind in 1:n_lambdas){
      for (l2_ind in 1:n_lambdas){
        fit <- try(
          learner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen, 
                  lambda_1 = lambda_1_vals[l1_ind], 
                  lambda_2 = lambda_2_vals[l2_ind], 
                  step_size = c, 
                  control = list(max_iter = max_iter)))
        if ('try-error' %in% class(fit)){
          norms_temp <- NA
        } else {
          norms_temp <- sum((fit$learner_estimate - theta_hat_0_external)^2)
        }
        norms_training[l1_ind, l2_ind] <- norms_temp
      }
    }
    min_ind <- which(norms_training == min(norms_training, na.rm = TRUE), arr.ind = TRUE)
    fit <- learner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen, 
                   lambda_1 = lambda_1_vals[min_ind[1]], 
                   lambda_2 = lambda_2_vals[min_ind[2]], 
                   step_size = c, 
                   control = list(max_iter = max_iter))
    res_learner <- norm(fit$learner_estimate - thetas$theta_0, type = 'F')
    
    ## Plotting LEARNER objective function
    pdf(paste0('ObjFunc-', title_help, '_var_1=', round(var_1, 2),'.pdf'), width = 6, height = 4)
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    plot(fit$objective_values, type = 'l', col = 'blue', lwd = 2, 
         xlab = 'Iteration number', ylab = 'Objective function')
    dev.off()
    
    return(list(tsvd = res_tsvd,
                dlearner = res_dlearner, learner = res_learner,
                lambda_1_ind = min_ind[1], 
                lambda_2_ind = min_ind[2], 
                r_chosen = r_chosen))
  }
  tsvd <- mean(sapply(res, '[[', 'tsvd'))
  dlearner <- mean(sapply(res, '[[', 'dlearner'))
  learner <- mean(sapply(res, '[[', 'learner'))
  
  tsvd_sd <- sd(sapply(res, '[[', 'tsvd'))
  dlearner_sd <- sd(sapply(res, '[[', 'dlearner'))
  learner_sd <- sd(sapply(res, '[[', 'learner'))
  
  lambda_best <- matrix(0, nrow = n_lambdas, ncol = n_lambdas)
  for (i in 1:n_iter){
    lambda_best[res[[i]]$lambda_1_ind, res[[i]]$lambda_2_ind] <-
      lambda_best[res[[i]]$lambda_1_ind, res[[i]]$lambda_2_ind] + 1
  }
  lambda_best <- lambda_best / n_iter
  r_chosen <- rep(0, times = max_rank)
  for (i in 1:n_iter){
    r_chosen[res[[i]]$r_chosen] <- r_chosen[res[[i]]$r_chosen] + 1
  }
  
  return(list(
    tsvd = tsvd, dlearner = dlearner, learner = learner,
    tsvd_sd = tsvd_sd, dlearner_sd = dlearner_sd, learner_sd = learner_sd,
    lambda_best = lambda_best,  r_chosen = r_chosen
  ))
}