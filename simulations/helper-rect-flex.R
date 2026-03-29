################################################################################
## learnerv2 (row/column penalties)
################################################################################

ensure_learnerv2_v03 <- function(tar_dir = getwd()){
  tar_v03 <- file.path(tar_dir, 'learnerv2_0.3.0.tar.gz')
  if (!file.exists(tar_v03)){
    stop('Missing learnerv2 tarball; run from simulations/ or pass tar_dir=')
  }
  if (!('learnerv2' %in% rownames(installed.packages())) || packageVersion('learnerv2') < '0.3.0'){
    message('Installing learnerv2 from: ', tar_v03)
    install.packages(tar_v03, repos = NULL, type = 'source')
  }
  if ('package:learner' %in% search()){
    detach('package:learner', unload = TRUE, character.only = TRUE)
  }
  library('learnerv2')
  invisible(TRUE)
}

################################################################################
## Simulation driver
## Relies on globals (set in runsim before calling), like helper.R::run_all_scenarios:
##   lambda_1_vals, lambda_2_vals, c (step size), max_iter, n_iter, max_rank, title_help
## Also uses helper.R globals: p, q, set_theta, get_theta_hat
################################################################################

# TSVD, D-LEARNER, LEARNER (tied), LEARNER (separate row/col penalties)
run_all_scenarios_rect_flex <- function(var_0_all, var_1_all, r, same_latent_space = TRUE, dif_factor = 4,
                                        independent_noise = TRUE){
  ensure_learnerv2_v03()

  n_lambdas <- length(lambda_1_vals)
  if (length(lambda_2_vals) != n_lambdas){
    stop('lambda_1_vals and lambda_2_vals must have the same length')
  }

  max_rank_int <- max(as.integer(floor(max_rank)), 1L)

  thetas <- set_theta(r = r, same_latent_space = same_latent_space, dif_factor = dif_factor)

  run_1_scenario <- function(var_0, var_1){
    res <- foreach(i = 1:n_iter, .packages = c('learnerv2', 'ScreeNOT')) %dorng% {
      ## Data Generation
      theta_hat_0 <- get_theta_hat(thetas$theta_0, var_val = var_0, independent_noise = independent_noise)
      theta_hat_1 <- get_theta_hat(thetas$theta_1, var_val = var_1, independent_noise = independent_noise)

      ## Truncated SVD
      svd0 <- svd(theta_hat_0, nu = r, nv = r)
      tsvd_est <- svd0$u %*% diag(svd0$d[1:r], nrow = r, ncol = r) %*% t(svd0$v)
      res_tsvd = norm(tsvd_est - thetas$theta_0, type = 'F')

      ## D-LEARNER
      r_chosen <- max(adaptiveHardThresholding(Y = theta_hat_1, k = max_rank_int)$r, 1)
      dlearner_est <- dlearner(Y_source = theta_hat_1, Y_target = theta_hat_0, r = r_chosen)$dlearner_estimate
      res_dlearner = norm(dlearner_est - thetas$theta_0, type = 'F')

      ## LEARNER (tied row/col penalties)
      cv_tied <- cv.learner(Y_source = theta_hat_1, Y_target = theta_hat_0,
                            r_source = r_chosen, r_target = r_chosen,
                            lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals,
                            step_size = c,
                            control = list(max_iter = max_iter))
      fit_tied <- learner(Y_source = theta_hat_1, Y_target = theta_hat_0,
                          r_source = r_chosen, r_target = r_chosen,
                          lambda_1 = cv_tied$lambda_1_min, lambda_2 = cv_tied$lambda_2_min,
                          step_size = c,
                          control = list(max_iter = max_iter))
      res_learner_tied <- norm(fit_tied$learner_estimate - thetas$theta_0, type = 'F')

      lambda_1_tied_ind <- which(lambda_1_vals == cv_tied$lambda_1_min)[1]
      lambda_2_tied_ind <- which(lambda_2_vals == cv_tied$lambda_2_min)[1]

      ## Plotting LEARNER objective function (tied penalties)
      pdf(paste0('ObjFunc-', title_help, '_tied_var_1=', round(var_1, 2), '_rep=', i, '.pdf'), width = 6, height = 4)
      par(mar = c(5.1, 4.6, 4.1, 2.1))
      plot(fit_tied$objective_values, type = 'l', col = 'blue', lwd = 2,
           xlab = 'Iteration number', ylab = 'Objective function')
      dev.off()

      ## LEARNER (separate row/col penalties)
      cv_sep <- cv.learner(Y_source = theta_hat_1, Y_target = theta_hat_0,
                           r_source = r_chosen, r_target = r_chosen,
                           lambda_1_row_all = lambda_1_vals, lambda_1_col_all = lambda_1_vals,
                           lambda_2_all = lambda_2_vals,
                           step_size = c,
                           control = list(max_iter = max_iter))
      fit_sep <- learner(Y_source = theta_hat_1, Y_target = theta_hat_0,
                         r_source = r_chosen, r_target = r_chosen,
                         lambda_1_row = cv_sep$lambda_1_row_min, lambda_1_col = cv_sep$lambda_1_col_min,
                         lambda_2 = cv_sep$lambda_2_min,
                         step_size = c,
                         control = list(max_iter = max_iter))
      res_learner_sep <- norm(fit_sep$learner_estimate - thetas$theta_0, type = 'F')

      lambda_1_row_ind <- which(lambda_1_vals == cv_sep$lambda_1_row_min)[1]
      lambda_1_col_ind <- which(lambda_1_vals == cv_sep$lambda_1_col_min)[1]
      lambda_2_sep_ind <- which(lambda_2_vals == cv_sep$lambda_2_min)[1]

      ## Plotting LEARNER objective function (separate row/col penalties)
      pdf(paste0('ObjFunc-', title_help, '_sep_var_1=', round(var_1, 2), '_rep=', i, '.pdf'), width = 6, height = 4)
      par(mar = c(5.1, 4.6, 4.1, 2.1))
      plot(fit_sep$objective_values, type = 'l', col = 'blue', lwd = 2,
           xlab = 'Iteration number', ylab = 'Objective function')
      dev.off()

      return(list(tsvd = res_tsvd,
                  dlearner = res_dlearner,
                  learner_tied = res_learner_tied,
                  learner_sep = res_learner_sep,
                  lambda_1_tied_ind = lambda_1_tied_ind,
                  lambda_2_tied_ind = lambda_2_tied_ind,
                  lambda_1_row_ind = lambda_1_row_ind,
                  lambda_1_col_ind = lambda_1_col_ind,
                  lambda_2_sep_ind = lambda_2_sep_ind,
                  r_chosen = r_chosen))
    }

    tsvd <- mean(sapply(res, '[[', 'tsvd'))
    tsvd_sd <- sd(sapply(res, '[[', 'tsvd'))
    dlearner <- mean(sapply(res, '[[', 'dlearner'))
    dlearner_sd <- sd(sapply(res, '[[', 'dlearner'))
    learner_tied <- mean(sapply(res, '[[', 'learner_tied'))
    learner_tied_sd <- sd(sapply(res, '[[', 'learner_tied'))
    learner_sep <- mean(sapply(res, '[[', 'learner_sep'))
    learner_sep_sd <- sd(sapply(res, '[[', 'learner_sep'))
    r_chosen_mean <- mean(sapply(res, '[[', 'r_chosen'))

    lambda_best_tied <- matrix(0, nrow = n_lambdas, ncol = n_lambdas)
    for (k in 1:n_iter){
      lambda_best_tied[res[[k]]$lambda_1_tied_ind, res[[k]]$lambda_2_tied_ind] <-
        lambda_best_tied[res[[k]]$lambda_1_tied_ind, res[[k]]$lambda_2_tied_ind] + 1
    }
    lambda_best_tied <- lambda_best_tied / n_iter

    lambda_best_sep <- array(0, dim = c(n_lambdas, n_lambdas, n_lambdas))
    for (k in 1:n_iter){
      lambda_best_sep[res[[k]]$lambda_1_row_ind, res[[k]]$lambda_1_col_ind, res[[k]]$lambda_2_sep_ind] <-
        lambda_best_sep[res[[k]]$lambda_1_row_ind, res[[k]]$lambda_1_col_ind, res[[k]]$lambda_2_sep_ind] + 1
    }
    lambda_best_sep <- lambda_best_sep / n_iter

    r_chosen_hist <- rep(0, times = max_rank_int)
    for (k in 1:n_iter){
      rc <- res[[k]]$r_chosen
      if (rc <= length(r_chosen_hist)){
        r_chosen_hist[rc] <- r_chosen_hist[rc] + 1
      }
    }

    return(list(tsvd = tsvd, tsvd_sd = tsvd_sd,
                dlearner = dlearner, dlearner_sd = dlearner_sd,
                learner_tied = learner_tied, learner_tied_sd = learner_tied_sd,
                learner_sep = learner_sep, learner_sep_sd = learner_sep_sd,
                r_chosen_mean = r_chosen_mean,
                lambda_best_tied = lambda_best_tied,
                lambda_best_sep = lambda_best_sep,
                r_chosen = r_chosen_hist))
  }

  n_scenarios <- length(var_1_all)
  res <- vector(mode = 'list', length = n_scenarios)

  for (j in 1:n_scenarios){
    print(paste('Running var_1=', var_1_all[j]))
    res[[j]] <- run_1_scenario(var_0 = var_0_all, var_1 = var_1_all[j])
  }

  return(list(tsvd = sapply(res, '[[', 'tsvd'),
              tsvd_sd = sapply(res, '[[', 'tsvd_sd'),
              dlearner = sapply(res, '[[', 'dlearner'),
              dlearner_sd = sapply(res, '[[', 'dlearner_sd'),
              learner_tied = sapply(res, '[[', 'learner_tied'),
              learner_tied_sd = sapply(res, '[[', 'learner_tied_sd'),
              learner_sep = sapply(res, '[[', 'learner_sep'),
              learner_sep_sd = sapply(res, '[[', 'learner_sep_sd'),
              r_chosen_mean = sapply(res, '[[', 'r_chosen_mean'),
              lambda_best_tied = lapply(res, '[[', 'lambda_best_tied'),
              lambda_best_sep = lapply(res, '[[', 'lambda_best_sep'),
              r_chosen = lapply(res, '[[', 'r_chosen'),
              p = p, q = q, r_true = r,
              same_latent_space = same_latent_space,
              dif_factor = dif_factor,
              var_0_all = var_0_all, var_1_all = var_1_all,
              n_iter = n_iter, max_rank = max_rank_int,
              c = c,
              lambda_1_all = lambda_1_vals, lambda_2_all = lambda_2_vals))
}
