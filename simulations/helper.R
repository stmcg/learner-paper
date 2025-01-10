################################################################################
## Loading packages and libraries
################################################################################

library('MASS')
library('doParallel')
library('doRNG')
library('foreach')
library('ScreeNOT')
library('learner')

################################################################################
## Global parameters
################################################################################

n_cores <- 10
registerDoParallel(cores=n_cores)

q <- 50
p <- 5000
n_iter <- 50
var_0_all <- 0.1
var_1_all <- var_0_all / c(1, 3, 5, 10)
n_scenarios <- length(var_1_all)
max_rank <- min(q, p) / 3

################################################################################
## Helper functions
################################################################################

# Obtains theta_0 and theta_1
set_theta <- function(r, same_latent_space, dif_factor){
  theta_temp <- matrix(rnorm(n = p * q, mean = 0, sd = 1), nrow = p, ncol = q)
  theta_svd <- svd(theta_temp, nu = r, nv = r)
  sing_vals <- theta_svd$d[1:r]
  theta_0 <- theta_svd$u %*% diag(sing_vals, nrow = r, ncol = r) %*% t(theta_svd$v)

  if (same_latent_space){
    u_new <- theta_svd$u
    v_new <- theta_svd$v
  } else {
    threshold_p <- 1/(dif_factor * sqrt(p)); threshold_q <- 1/(dif_factor * sqrt(q))
    delta1 <- matrix(runif(p * r, min = -threshold_p, max = threshold_p), nrow = p, ncol = r)
    delta2 <- matrix(runif(q * r, min = -threshold_q, max = threshold_q), nrow = q, ncol = r)
    u_new <-  qr.Q(qr(theta_svd$u + delta1))
    v_new <- qr.Q(qr(theta_svd$v + delta2))
  }
  theta_1 <- u_new %*% diag(rev(sing_vals), nrow = r, ncol = r) %*% t(v_new)

  print(paste('d_U = ', norm(u_new %*% t(u_new) - theta_svd$u %*% t(theta_svd$u), type = 'F')))
  print(paste('d_V = ', norm(v_new %*% t(v_new) - theta_svd$v %*% t(theta_svd$v), type = 'F')))
  return(list(theta_0 = theta_0, theta_1 = theta_1))
}

# Obtains Y_0 and Y_1
get_theta_hat <- function(theta, var_val, independent_noise){
  if (independent_noise){
    theta_hat <- theta +
      matrix(rnorm(p * q, mean = 0, sd = sqrt(var_val)), nrow = p, ncol = q)
  } else {
    Sigma <- matrix(rho, nrow = p, ncol = p)
    diag(Sigma) <- 1
    Sigma <- Sigma * var_val
    theta_hat <- theta + 
      t(mvrnorm(n = q, mu = rep(0, length = p), Sigma = Sigma))
  }
  return(theta_hat)
}

# Performs the simulation study for one data generation scenario
run_1_scenario <- function(var_0, var_1, r, thetas, independent_noise){
  res <- foreach(i = 1:n_iter) %dorng% {
    ## Data Generation
    theta_hat_0 <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
    theta_hat_1 <- get_theta_hat(thetas$theta_1, var_val = var_1, independent_noise = independent_noise)

    ## Truncated SVD
    svd0 <- svd(theta_hat_0, nu = r, nv = r)
    tsvd_est <- svd0$u %*% diag(svd0$d[1:r], nrow = r, ncol = r) %*% t(svd0$v)
    res_tsvd = norm(tsvd_est - thetas$theta_0, type = 'F')
    
    ## D-LEARNER
    r_chosen <- max(adaptiveHardThresholding(Y = theta_hat_1, k = max_rank)$r, 1)
    dlearner_est <- dlearner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen)$dlearner_estimate
    res_dlearner = norm(dlearner_est - thetas$theta_0, type = 'F')

    ## LEARNER
    res_cv <- cv.learner(Y_source = theta_hat_1, Y_target = theta_hat_0, 
                         r = r_chosen, 
                         lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals, 
                         step_size = c, 
                         control = list(max_iter = max_iter))
    lambda_1_min <- res_cv$lambda_1_min; lambda_2_min <- res_cv$lambda_2_min
    fit <- learner(Y_source = theta_hat_1, Y_target = theta_hat_0, 
                   r = r_chosen, lambda_1 = lambda_1_min, lambda_2 = lambda_2_min, 
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
                lambda_1_ind = which(lambda_1_vals == lambda_1_min),
                lambda_2_ind = which(lambda_2_vals == lambda_2_min),
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

# Performs the simulation study for all data generation scenarios for a fixed rank
run_all_scenarios <- function(var_0_all, var_1_all, r, same_latent_space = TRUE, dif_factor = 4, 
                              independent_noise = TRUE){
  theta_temp <- set_theta(r = r, same_latent_space = same_latent_space, dif_factor = dif_factor)
  res_tsvd <- res_dlearner <- res_learner <-
    res_tsvd_sd <- res_dlearner_sd <- res_learner_sd <- rep(NA, times = n_scenarios)
  lambda_best <- r_chosen <- vector(mode = 'list', length = n_scenarios)
  for (i in 1:n_scenarios){
    temp <- run_1_scenario(var_0 = var_0_all, var_1 = var_1_all[i], r = r, thetas = theta_temp, 
                           independent_noise = independent_noise)
    res_tsvd[i] <- temp$tsvd
    res_dlearner[i] <- temp$dlearner
    res_learner[i] <- temp$learner

    res_tsvd_sd[i] <- temp$tsvd_sd
    res_dlearner_sd[i] <- temp$dlearner_sd
    res_learner_sd[i] <- temp$learner_sd

    lambda_best[[i]] <- temp$lambda_best
    r_chosen[[i]] <- temp$r_chosen
  }
  return(list(tsvd = res_tsvd, dlearner = res_dlearner, learner = res_learner,
              tsvd_sd = res_tsvd_sd, dlearner_sd = res_dlearner_sd, learner_sd = res_learner_sd,
              lambda_best = lambda_best, r_chosen = r_chosen))
}
