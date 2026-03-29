################################################################################
## Multisource: Source A high similarity to target, Source B stronger low similarity
## (same DGP as runsim-same vs runsim-dif-verylow: dif_factor = 2 for source B perturbation).
## Equal source noise variance (0.01 each). Target noise var_0 = 0.1 (primary).
## Independent noise; p x q = 5000 x 50; r = 4.
##
## Tuning:
##   Source A only: lambda / step_size as runsim-same.R (high-sim primary).
##   Source B only: lambda / step_size as runsim-dif-verylow.R (matches DGP dif_factor = 2).
##   List (two sources): wider lambda grid; step size c_list = midpoint of c_same and c_dif_verylow.
##
## Weights for weighted list/D-LEARNER: inverse squared subspace distance to target
## (helper-multisource-high-low.R).
##
## learnerv2 from local learnerv2*.tar.gz (same as runsim-multisource.R).
################################################################################

rm(list = ls())

if (!('learnerv2' %in% rownames(installed.packages()))){
  tar <- list.files(getwd(), pattern = '^learnerv2.*\\.tar\\.gz$', full.names = TRUE)
  if (length(tar) == 0) stop('learnerv2 not installed and no learnerv2*.tar.gz found in current directory')
  install.packages(tar[1], repos = NULL, type = 'source')
}

source('helper.R')

detach('package:learner', unload = TRUE, character.only = TRUE)
library('learnerv2')
source('helper-multisource-high-low.R')

title_help <- 'multisource-high-low'

## Primary grids (must match lengths for cv.learner 2D grids)
n_lambdas <- 5

## runsim-same.R — source A (high similarity)
lambda_1_vals_a <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
lambda_2_vals_a <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c_a <- 0.0035

## runsim-dif-verylow.R — source B (low similarity; same lambdas as dif-moderate, c = 0.07)
lambda_1_vals_b <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals_b <- 10^seq(from = -2, to = 1, length.out = n_lambdas)
c_b <- 0.07

## List: wide lambda grid (more points, full range covering both single-source regimes)
lambda_1_vals_list <- 10^seq(from = 0, to = 4, length.out = n_lambdas)
lambda_2_vals_list <- 10^seq(from = -2, to = 1, length.out = n_lambdas)
c_list <- (c_a + c_b) / 2

max_iter <- 75

var_0 <- 0.1
var_1_a <- 0.01
var_1_b <- 0.01

r <- 4
dif_factor <- 2

set.seed(12345)
thetas <- set_theta_high_low_sources(r = r, dif_factor = dif_factor)

set.seed(12345)
res_multisource_high_low <- run_1_scenario_two_sources_high_low(
  var_0 = var_0,
  var_1_a = var_1_a,
  var_1_b = var_1_b,
  r = r,
  thetas = thetas,
  independent_noise = TRUE
)

save(res_multisource_high_low, thetas,
     file = 'simres-multisource-high-low.RData', compress = 'xz')
