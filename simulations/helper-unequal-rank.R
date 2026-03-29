################################################################################
## Helper functions for unequal-rank simulations
################################################################################

# Obtains theta_0 and theta_1 when the target and source have different ranks
set_theta_unequal_rank <- function(r_0, r_1){
  r_max <- max(r_0, r_1)
  
  # Generate base SVD with enough vectors for both ranks
  theta_temp <- matrix(rnorm(n = p * q, mean = 0, sd = 1), nrow = p, ncol = q)
  theta_svd <- svd(theta_temp, nu = r_max, nv = r_max)
  
  # theta_0 uses first r_0 singular vectors/values
  u_0 <- theta_svd$u[, 1:r_0, drop = FALSE]
  v_0 <- theta_svd$v[, 1:r_0, drop = FALSE]
  sing_vals_0 <- theta_svd$d[1:r_0]
  theta_0 <- u_0 %*% diag(sing_vals_0, nrow = r_0, ncol = r_0) %*% t(v_0)
  
  # theta_1 uses first r_1 singular vectors/values
  # When r_1 > r_0: shares u_0/v_0 exactly, plus extra orthogonal directions
  # When r_1 < r_0: uses only a subset of u_0/v_0
  r_shared <- min(r_0, r_1)
  u_1 <- theta_svd$u[, 1:r_1, drop = FALSE]
  v_1 <- theta_svd$v[, 1:r_1, drop = FALSE]
  sing_vals_1 <- theta_svd$d[1:r_1]
  sing_vals_1_rev <- c(rev(sing_vals_1[1:r_shared]), sing_vals_1[(r_shared+1):r_1])
  theta_1 <- u_1 %*% diag(sing_vals_1_rev, nrow = r_1, ncol = r_1) %*% t(v_1)
  
  # Diagnostic: similarity over the shared directions
  pu_0 <- u_0[, 1:r_shared, drop = FALSE] %*% t(u_0[, 1:r_shared, drop = FALSE])
  pu_1 <- u_1[, 1:r_shared, drop = FALSE] %*% t(u_1[, 1:r_shared, drop = FALSE])
  pv_0 <- v_0[, 1:r_shared, drop = FALSE] %*% t(v_0[, 1:r_shared, drop = FALSE])
  pv_1 <- v_1[, 1:r_shared, drop = FALSE] %*% t(v_1[, 1:r_shared, drop = FALSE])
  print(paste('d_U =', norm(pu_1 - pu_0, type = 'F')))
  print(paste('d_V =', norm(pv_1 - pv_0, type = 'F')))
  
  return(list(theta_0 = theta_0, theta_1 = theta_1))
}

# Performs the simulation study for one data generation scenario
run_1_scenario_unequal_rank <- function(var_0, var_1, r_all, r_0, thetas, independent_noise){
  max_r <- max(r_all)
  res <- foreach(i = 1:n_iter) %dorng% {
    ## Data Generation
    theta_hat_0 <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
    theta_hat_1 <- get_theta_hat(thetas$theta_1, var_val = var_1, independent_noise = independent_noise)
    
    ## Truncated SVD (target population): apply once before loop, using correct target rank r_0
    svd0 <- svd(theta_hat_0, nu = max_r, nv = max_r)
    tsvd_est_0 <- svd0$u[, 1:r_0, drop = FALSE] %*%
      diag(svd0$d[1:r_0], nrow = r_0, ncol = r_0) %*%
      t(svd0$v[, 1:r_0, drop = FALSE])
    res_tsvd_0 <- norm(tsvd_est_0 - thetas$theta_0, type = 'F')
    
    out <- vector(mode = 'list', length = length(r_all))
    names(out) <- paste0('r', r_all)
    if (r_0 == 4){
      c <- c_r0_4
    } else {
      c <- c_r0_8
    }
    for (r_ind in seq_along(r_all)){
      r <- r_all[r_ind]
      ## D-LEARNER
      dlearner_est <- dlearner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r)$dlearner_estimate
      res_dlearner = norm(dlearner_est - thetas$theta_0, type = 'F')
      
      ## LEARNER
      res_cv <- cv.learner(Y_source = theta_hat_1, Y_target = theta_hat_0, 
                           r = r, 
                           lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals, 
                           step_size = c, 
                           control = list(max_iter = max_iter))
      lambda_1_min <- res_cv$lambda_1_min; lambda_2_min <- res_cv$lambda_2_min
      fit <- learner(Y_source = theta_hat_1, Y_target = theta_hat_0, 
                     r = r, lambda_1 = lambda_1_min, lambda_2 = lambda_2_min, 
                     step_size = c, 
                     control = list(max_iter = max_iter))
      res_learner <- norm(fit$learner_estimate - thetas$theta_0, type = 'F')
      
      ## Plotting LEARNER objective function
      pdf(paste0('ObjFunc-', title_help, '_r=', r, '_var_1=', round(var_1, 2),'.pdf'), width = 6, height = 4)
      par(mar = c(5.1, 4.6, 4.1, 2.1))
      plot(fit$objective_values, type = 'l', col = 'blue', lwd = 2, 
           xlab = 'Iteration number', ylab = 'Objective function')
      dev.off()
      
      out[[r_ind]] <- list(tsvd = NA_real_,
                           dlearner = res_dlearner, learner = res_learner,
                           lambda_1_ind = which(lambda_1_vals == lambda_1_min),
                           lambda_2_ind = which(lambda_2_vals == lambda_2_min))
    }
    ## Fill TSVD into the slot for target rank r_0
    r_0_ind <- which(r_all == r_0)
    if (length(r_0_ind) > 0) out[[r_0_ind]]$tsvd <- res_tsvd_0
    return(out)
  }
  out_all <- vector(mode = 'list', length = length(r_all))
  names(out_all) <- paste0('r', r_all)
  for (r_ind in seq_along(r_all)){
    res_rank <- lapply(res, '[[', r_ind)
    tsvd <- mean(sapply(res_rank, '[[', 'tsvd'), na.rm = TRUE)
    dlearner <- mean(sapply(res_rank, '[[', 'dlearner'))
    learner <- mean(sapply(res_rank, '[[', 'learner'))
    
    tsvd_sd <- sd(sapply(res_rank, '[[', 'tsvd'), na.rm = TRUE)
    dlearner_sd <- sd(sapply(res_rank, '[[', 'dlearner'))
    learner_sd <- sd(sapply(res_rank, '[[', 'learner'))
    
    lambda_best <- matrix(0, nrow = n_lambdas, ncol = n_lambdas)
    for (i in 1:n_iter){
      lambda_best[res_rank[[i]]$lambda_1_ind, res_rank[[i]]$lambda_2_ind] <-
        lambda_best[res_rank[[i]]$lambda_1_ind, res_rank[[i]]$lambda_2_ind] + 1
    }
    lambda_best <- lambda_best / n_iter
    
    out_all[[r_ind]] <- list(
      tsvd = tsvd, dlearner = dlearner, learner = learner,
      tsvd_sd = tsvd_sd, dlearner_sd = dlearner_sd, learner_sd = learner_sd,
      lambda_best = lambda_best
    )
  }
  return(out_all)
}

# Performs the simulation study for all data generation scenarios for fixed ranks
run_all_scenarios_unequal_rank <- function(var_0_all, var_1_all, r_0, r_1, 
                                           r_all = c(4, 8), independent_noise = TRUE){
  theta_temp <- set_theta_unequal_rank(r_0 = r_0, r_1 = r_1)
  out_all <- vector(mode = 'list', length = length(r_all))
  names(out_all) <- paste0('r', r_all)
  for (r_ind in seq_along(r_all)){
    out_all[[r_ind]] <- list(
      tsvd = rep(NA, times = n_scenarios),
      dlearner = rep(NA, times = n_scenarios),
      learner = rep(NA, times = n_scenarios),
      tsvd_sd = rep(NA, times = n_scenarios),
      dlearner_sd = rep(NA, times = n_scenarios),
      learner_sd = rep(NA, times = n_scenarios),
      lambda_best = vector(mode = 'list', length = n_scenarios)
    )
  }
  for (i in 1:n_scenarios){
    temp <- run_1_scenario_unequal_rank(var_0 = var_0_all, var_1 = var_1_all[i], r_all = r_all, r_0 = r_0,
                                        thetas = theta_temp, independent_noise = independent_noise)
    for (r_ind in seq_along(r_all)){
      out_all[[r_ind]]$tsvd[i] <- temp[[r_ind]]$tsvd
      out_all[[r_ind]]$dlearner[i] <- temp[[r_ind]]$dlearner
      out_all[[r_ind]]$learner[i] <- temp[[r_ind]]$learner
      
      out_all[[r_ind]]$tsvd_sd[i] <- temp[[r_ind]]$tsvd_sd
      out_all[[r_ind]]$dlearner_sd[i] <- temp[[r_ind]]$dlearner_sd
      out_all[[r_ind]]$learner_sd[i] <- temp[[r_ind]]$learner_sd
      
      out_all[[r_ind]]$lambda_best[[i]] <- temp[[r_ind]]$lambda_best
    }
  }
  return(out_all)
}
