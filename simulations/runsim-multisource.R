################################################################################
## Independent noise, rank 4, high similarity.
## Two Y_source: first has source noise var = 0.1, second has var = 0.01 (1/10).
## Uses learnerv2 (multi-source) from local tar.gz in current working directory.
################################################################################

rm(list = ls())

## Use learnerv2 for multi-source LEARNER/D-LEARNER (install from tar.gz if needed)
## Expected file in cwd: learnerv2_0.2.0.tar.gz (or any learnerv2_*.tar.gz)
if (!("learnerv2" %in% rownames(installed.packages()))) {
  tar <- list.files(getwd(), pattern = "^learnerv2.*\\.tar\\.gz$", full.names = TRUE)
  if (length(tar) == 0) stop("learnerv2 not installed and no learnerv2*.tar.gz found in current directory")
  install.packages(tar[1], repos = NULL, type = "source")
}
source('helper.R')
detach("package:learner", unload = TRUE)
library(learnerv2)
source('helper-multisource.R')

title_help <- 'multisource'

# Target variance; two source variances (second = first / 10)
var_0 <- 0.1
var_1_a <- 0.01
var_1_b <- 0.5
## Weighted list/D-LEARNER: per-rep subspace weights (see helper-multisource.R; same rule as high-low).

# Same tuning as runsim-same.R
n_lambdas <- 5
lambda_1_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
lambda_2_vals <- 10^seq(from = -4, to = 4, length.out = n_lambdas)
c <- 0.0045
max_iter <- 75

set.seed(12345)
thetas <- set_theta(r = 4, same_latent_space = TRUE, dif_factor = 4)
res_indep_r4_high_0.1 <- run_1_scenario_two_sources(
  var_0 = var_0,
  var_1_a = var_1_a,
  var_1_b = var_1_b,
  r = 4,
  thetas = thetas,
  independent_noise = TRUE
)

save.image(file = 'simres-multisource.RData')
