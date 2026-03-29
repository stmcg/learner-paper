################################################################################
## One scenario with two Y_source: same theta_1, different noise variances.
## Second source has variance = (first source variance) / 10.
## Uses get_theta_hat, lambda_*_vals, c, max_iter, n_iter, max_rank from helper.R
##
## Weighted multisource (list-w): same construction as helper-multisource-high-low.R:
## d = ||P_Uk-P_U0||_F^2+||P_Vk-P_V0||_F^2 with rank-r SVDs of Y_0 and each source;
## w_k ∝ 1/(d_k+eps), normalized (per replication).
################################################################################

# d_{U,k} + d_{V,k} = ||P_Uk - P_U0||_F^2 + ||P_Vk - P_V0||_F^2 with rank-r_trunc SVDs
subspace_dist_sq_uv <- function(Y0, Yk, r_trunc){
  s0 <- svd(Y0, nu = r_trunc, nv = r_trunc)
  sk <- svd(Yk, nu = r_trunc, nv = r_trunc)
  U0 <- s0$u[, 1:r_trunc, drop = FALSE]
  V0 <- s0$v[, 1:r_trunc, drop = FALSE]
  Uk <- sk$u[, 1:r_trunc, drop = FALSE]
  Vk <- sk$v[, 1:r_trunc, drop = FALSE]
  PU0 <- U0 %*% t(U0)
  PUk <- Uk %*% t(Uk)
  PV0 <- V0 %*% t(V0)
  PVk <- Vk %*% t(Vk)
  norm(PUk - PU0, type = 'F')^2 + norm(PVk - PV0, type = 'F')^2
}

weights_from_subspace_d <- function(d_a, d_b, eps = 1e-12){
  w_raw <- c(1 / (d_a + eps), 1 / (d_b + eps))
  w_raw / sum(w_raw)
}

run_1_scenario_two_sources <- function(var_0, var_1_a, var_1_b, r, thetas, independent_noise) {
  res <- foreach(i = 1:n_iter) %dorng% {
    ## Data generation: one target, two sources (different noise levels)
    theta_hat_0 <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
    theta_hat_1a <- get_theta_hat(thetas$theta_1, var_val = var_1_a, independent_noise = independent_noise)
    theta_hat_1b <- get_theta_hat(thetas$theta_1, var_val = var_1_b, independent_noise = independent_noise)

    ## Truncated SVD (target only; same for both sources)
    svd0 <- svd(theta_hat_0, nu = r, nv = r)
    tsvd_est <- svd0$u %*% diag(svd0$d[1:r], nrow = r, ncol = r) %*% t(svd0$v)
    res_tsvd <- norm(tsvd_est - thetas$theta_0, type = 'F')

    ## D-LEARNER with source A
    r_chosen_a <- max(adaptiveHardThresholding(Y = theta_hat_1a, k = max_rank)$r, 1)
    dlearner_est_a <- dlearner(Y_source = theta_hat_1a, Y_target = theta_hat_0, r = r_chosen_a)$dlearner_estimate
    res_dlearner_a <- norm(dlearner_est_a - thetas$theta_0, type = 'F')

    ## D-LEARNER with source B
    r_chosen_b <- max(adaptiveHardThresholding(Y = theta_hat_1b, k = max_rank)$r, 1)
    dlearner_est_b <- dlearner(Y_source = theta_hat_1b, Y_target = theta_hat_0, r = r_chosen_b)$dlearner_estimate
    res_dlearner_b <- norm(dlearner_est_b - thetas$theta_0, type = 'F')

    ## LEARNER with source A
    res_cv_a <- cv.learner(Y_source = theta_hat_1a, Y_target = theta_hat_0,
                           r_source = r_chosen_a, r_target = r_chosen_a,
                           lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals,
                           step_size = 0.001,
                           control = list(max_iter = max_iter))
    fit_a <- learner(Y_source = theta_hat_1a, Y_target = theta_hat_0,
                     r_source = r_chosen_a, r_target = r_chosen_a,
                     lambda_1 = res_cv_a$lambda_1_min, lambda_2 = res_cv_a$lambda_2_min,
                     step_size = 0.001,
                     control = list(max_iter = max_iter))
    res_learner_a <- norm(fit_a$learner_estimate - thetas$theta_0, type = 'F')
    pdf(paste0('ObjFunc-', title_help, '_source-a.pdf'), width = 6, height = 4)
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    plot(fit_a$objective_values, type = 'l', col = 'blue', lwd = 2,
         xlab = 'Iteration number', ylab = 'Objective function')
    dev.off()

    ## LEARNER with source B
    res_cv_b <- cv.learner(Y_source = theta_hat_1b, Y_target = theta_hat_0,
                           r_source = r_chosen_b, r_target = r_chosen_b,
                           lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals,
                           step_size = 0.01,
                           control = list(max_iter = max_iter))
    fit_b <- learner(Y_source = theta_hat_1b, Y_target = theta_hat_0,
                     r_source = r_chosen_b, r_target = r_chosen_b,
                     lambda_1 = res_cv_b$lambda_1_min, lambda_2 = res_cv_b$lambda_2_min,
                     step_size = 0.01,
                     control = list(max_iter = max_iter))
    res_learner_b <- norm(fit_b$learner_estimate - thetas$theta_0, type = 'F')
    pdf(paste0('ObjFunc-', title_help, '_source-b.pdf'), width = 6, height = 4)
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    plot(fit_b$objective_values, type = 'l', col = 'blue', lwd = 2,
         xlab = 'Iteration number', ylab = 'Objective function')
    dev.off()

    ## Rank for multi-source: use first source
    r_chosen_list <- r_chosen_a

    Y_source_list <- list(theta_hat_1a, theta_hat_1b)
    d_sub_a <- subspace_dist_sq_uv(theta_hat_0, theta_hat_1a, r)
    d_sub_b <- subspace_dist_sq_uv(theta_hat_0, theta_hat_1b, r)
    source_weights <- weights_from_subspace_d(d_sub_a, d_sub_b)
    print(paste0('rep=', i, ' w_A=', source_weights[1], ' w_B=', source_weights[2],
                 ' d_sub_A=', d_sub_a, ' d_sub_B=', d_sub_b))

    ## D-LEARNER with Y_source = list(theta_hat_1a, theta_hat_1b)
    dlearner_est_list <- dlearner(Y_source = Y_source_list, Y_target = theta_hat_0, r = r_chosen_list)$dlearner_estimate
    res_dlearner_list <- norm(dlearner_est_list - thetas$theta_0, type = 'F')

    ## D-LEARNER with subspace-based source_weights (same rule as high-low multisource)
    dlearner_est_list_w <- dlearner(Y_source = Y_source_list, Y_target = theta_hat_0, r = r_chosen_list,
                                    source_weights = source_weights)$dlearner_estimate
    res_dlearner_list_w <- norm(dlearner_est_list_w - thetas$theta_0, type = 'F')

    ## LEARNER with Y_source = list(theta_hat_1a, theta_hat_1b)
    res_cv_list <- cv.learner(Y_source = Y_source_list, Y_target = theta_hat_0,
                              r_source = r_chosen_list, r_target = r_chosen_list,
                              lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals,
                              step_size = 0.0035,
                              control = list(max_iter = max_iter))
    fit_list <- learner(Y_source = Y_source_list, Y_target = theta_hat_0,
                        r_source = r_chosen_list, r_target = r_chosen_list,
                        lambda_1 = res_cv_list$lambda_1_min, lambda_2 = res_cv_list$lambda_2_min,
                        step_size = 0.0035,
                        control = list(max_iter = max_iter))
    res_learner_list <- norm(fit_list$learner_estimate - thetas$theta_0, type = 'F')
    pdf(paste0('ObjFunc-', title_help, '_list.pdf'), width = 6, height = 4)
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    plot(fit_list$objective_values, type = 'l', col = 'blue', lwd = 2,
         xlab = 'Iteration number', ylab = 'Objective function')
    dev.off()

    ## LEARNER with list + subspace-based source_weights
    res_cv_list_w <- cv.learner(Y_source = Y_source_list, Y_target = theta_hat_0,
                                r_source = r_chosen_list, r_target = r_chosen_list,
                                lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals,
                                step_size = 0.001,
                                source_weights = source_weights,
                                control = list(max_iter = max_iter))
    fit_list_w <- learner(Y_source = Y_source_list, Y_target = theta_hat_0,
                          r_source = r_chosen_list, r_target = r_chosen_list,
                          lambda_1 = res_cv_list_w$lambda_1_min, lambda_2 = res_cv_list_w$lambda_2_min,
                          step_size = 0.001,
                          source_weights = source_weights,
                          control = list(max_iter = max_iter))
    res_learner_list_w <- norm(fit_list_w$learner_estimate - thetas$theta_0, type = 'F')
    pdf(paste0('ObjFunc-', title_help, '_list-w.pdf'), width = 6, height = 4)
    par(mar = c(5.1, 4.6, 4.1, 2.1))
    plot(fit_list_w$objective_values, type = 'l', col = 'blue', lwd = 2,
         xlab = 'Iteration number', ylab = 'Objective function')
    dev.off()

    return(list(
      tsvd = res_tsvd,
      dlearner_a = res_dlearner_a, dlearner_b = res_dlearner_b,
      dlearner_list = res_dlearner_list, dlearner_list_w = res_dlearner_list_w,
      learner_a = res_learner_a, learner_b = res_learner_b,
      learner_list = res_learner_list, learner_list_w = res_learner_list_w,
      r_chosen_a = r_chosen_a, r_chosen_b = r_chosen_b,
      d_sub_a = d_sub_a, d_sub_b = d_sub_b,
      w_a = source_weights[1], w_b = source_weights[2],
      lambda_1_ind_a = which(lambda_1_vals == res_cv_a$lambda_1_min),
      lambda_2_ind_a = which(lambda_2_vals == res_cv_a$lambda_2_min),
      lambda_1_ind_b = which(lambda_1_vals == res_cv_b$lambda_1_min),
      lambda_2_ind_b = which(lambda_2_vals == res_cv_b$lambda_2_min),
      lambda_1_ind_list = which(lambda_1_vals == res_cv_list$lambda_1_min),
      lambda_2_ind_list = which(lambda_2_vals == res_cv_list$lambda_2_min),
      lambda_1_ind_list_w = which(lambda_1_vals == res_cv_list_w$lambda_1_min),
      lambda_2_ind_list_w = which(lambda_2_vals == res_cv_list_w$lambda_2_min)
    ))
  }

  n_grid <- length(lambda_1_vals)
  if (length(lambda_2_vals) != n_grid) {
    stop('lambda_1_vals and lambda_2_vals must have the same length.')
  }

  lambda_best_a <- matrix(0, nrow = n_grid, ncol = n_grid)
  lambda_best_b <- matrix(0, nrow = n_grid, ncol = n_grid)
  lambda_best_list <- matrix(0, nrow = n_grid, ncol = n_grid)
  lambda_best_list_w <- matrix(0, nrow = n_grid, ncol = n_grid)
  for (j in 1:n_iter) {
    i1a <- res[[j]]$lambda_1_ind_a[1]
    i2a <- res[[j]]$lambda_2_ind_a[1]
    i1b <- res[[j]]$lambda_1_ind_b[1]
    i2b <- res[[j]]$lambda_2_ind_b[1]
    i1l <- res[[j]]$lambda_1_ind_list[1]
    i2l <- res[[j]]$lambda_2_ind_list[1]
    i1lw <- res[[j]]$lambda_1_ind_list_w[1]
    i2lw <- res[[j]]$lambda_2_ind_list_w[1]
    if (!is.na(i1a) && !is.na(i2a)) {
      lambda_best_a[i1a, i2a] <- lambda_best_a[i1a, i2a] + 1
    }
    if (!is.na(i1b) && !is.na(i2b)) {
      lambda_best_b[i1b, i2b] <- lambda_best_b[i1b, i2b] + 1
    }
    if (!is.na(i1l) && !is.na(i2l)) {
      lambda_best_list[i1l, i2l] <- lambda_best_list[i1l, i2l] + 1
    }
    if (!is.na(i1lw) && !is.na(i2lw)) {
      lambda_best_list_w[i1lw, i2lw] <- lambda_best_list_w[i1lw, i2lw] + 1
    }
  }
  lambda_best_a <- lambda_best_a / n_iter
  lambda_best_b <- lambda_best_b / n_iter
  lambda_best_list <- lambda_best_list / n_iter
  lambda_best_list_w <- lambda_best_list_w / n_iter

  list(
    tsvd = mean(sapply(res, '[[', 'tsvd')),
    tsvd_sd = sd(sapply(res, '[[', 'tsvd')),
    dlearner_a = mean(sapply(res, '[[', 'dlearner_a')),
    dlearner_a_sd = sd(sapply(res, '[[', 'dlearner_a')),
    dlearner_b = mean(sapply(res, '[[', 'dlearner_b')),
    dlearner_b_sd = sd(sapply(res, '[[', 'dlearner_b')),
    dlearner_list = mean(sapply(res, '[[', 'dlearner_list')),
    dlearner_list_sd = sd(sapply(res, '[[', 'dlearner_list')),
    dlearner_list_w = mean(sapply(res, '[[', 'dlearner_list_w')),
    dlearner_list_w_sd = sd(sapply(res, '[[', 'dlearner_list_w')),
    learner_a = mean(sapply(res, '[[', 'learner_a')),
    learner_a_sd = sd(sapply(res, '[[', 'learner_a')),
    learner_b = mean(sapply(res, '[[', 'learner_b')),
    learner_b_sd = sd(sapply(res, '[[', 'learner_b')),
    learner_list = mean(sapply(res, '[[', 'learner_list')),
    learner_list_sd = sd(sapply(res, '[[', 'learner_list')),
    learner_list_w = mean(sapply(res, '[[', 'learner_list_w')),
    learner_list_w_sd = sd(sapply(res, '[[', 'learner_list_w')),
    r_chosen_a = sapply(res, '[[', 'r_chosen_a'),
    r_chosen_b = sapply(res, '[[', 'r_chosen_b'),
    d_sub_a_mean = mean(sapply(res, '[[', 'd_sub_a')),
    d_sub_a_sd = sd(sapply(res, '[[', 'd_sub_a')),
    d_sub_b_mean = mean(sapply(res, '[[', 'd_sub_b')),
    d_sub_b_sd = sd(sapply(res, '[[', 'd_sub_b')),
    w_a_mean = mean(sapply(res, '[[', 'w_a')),
    w_a_sd = sd(sapply(res, '[[', 'w_a')),
    w_b_mean = mean(sapply(res, '[[', 'w_b')),
    w_b_sd = sd(sapply(res, '[[', 'w_b')),
    lambda_1_vals = lambda_1_vals,
    lambda_2_vals = lambda_2_vals,
    lambda_best_a = lambda_best_a,
    lambda_best_b = lambda_best_b,
    lambda_best_list = lambda_best_list,
    lambda_best_list_w = lambda_best_list_w,
    weighting_note = 'Weights = (1/(d+eps)) normalized; d = ||P_Uk-P_U0||_F^2+||P_Vk-P_V0||_F^2 with rank-r SVDs (helper-multisource-high-low.R)'
  )
}
