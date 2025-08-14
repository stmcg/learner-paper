################################################################################
## Loading data and libraries
################################################################################
rm(list = ls())
load('dat.RData')
load('apply_learner_dlearner/learner.RData')
load('apply_learner_dlearner/DLEARNER-estimate-BBJ.RData')
load('key.RData')

library('lattice')
library('RColorBrewer')

################################################################################
## Data processing
################################################################################

## Re-ordering columns of theta_hat matrices
theta_hat_0 <- theta_hat_0[, order(key$ICD_cat)]
theta_hat_1 <- theta_hat_1[, order(key$ICD_cat)]
fit$learner_estimate <- fit$learner_estimate[, order(key$ICD_cat)]
theta_hat_0_dlearner <- theta_hat_0_dlearner[, order(key$ICD_cat)]

## Computing SVD of theta_hat matrices
svd0_chosen <- svd(theta_hat_0, nu = r_chosen, nv = r_chosen)
svd1_chosen <- svd(theta_hat_1, nu = r_chosen, nv = r_chosen)
svd_learner <- svd(fit$learner_estimate, nu = r_chosen, nv = r_chosen)
svd_dlearner <- svd(theta_hat_0_dlearner, nu = r_chosen, nv = r_chosen)

## Computing projection matrices 
## Note: Memory issues can arise if PU_0, PU_1, PU_0_learner, and PU_0_dlearner
##       are all computed in the same session
PU_0 <- svd0_chosen$u %*% t(svd0_chosen$u)
PU_1 <- svd1_chosen$u %*% t(svd1_chosen$u)
PU_0_learner <- svd_learner$u %*% t(svd_learner$u)
PU_0_dlearner <- svd_dlearner$u %*% t(svd_dlearner$u)

set.seed(1234)
myind <- sort(sample.int(p, 500, replace = F))
PU_0_small <- PU_0[myind, myind]
PU_1_small <- PU_1[myind, myind]
PU_0_learner_small <- PU_0_learner[myind, myind]
PU_0_dlearner_small <- PU_0_dlearner[myind, myind]


################################################################################
## Plotting projection matrices onto the space of latent phenotypic factors
################################################################################

coul <- c(colorRampPalette(c("firebrick3", "white"))(10),
          colorRampPalette(c("white", "dodgerblue3"))(10))

myround <- function(x, min, max){
  return(pmax(pmin(x, max), min))
}

cex.main <- 2.25; cex.label <- 1.75; cex.scales <- 1.5
min_val <- -0.15; max_val <- 0.15

pdf('PV-source.pdf', width = 8, height = 8)
temp_source <- svd1_chosen$v %*% t(svd1_chosen$v)
levelplot(myround(temp_source, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Phenotype Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('Source Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales), y=list(cex=cex.scales)))
dev.off()

pdf('PV-target.pdf', width = 8, height = 8)
temp_target <- svd0_chosen$v %*% t(svd0_chosen$v)
levelplot(myround(temp_target, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Phenotype Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('Target Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

pdf('PV-target-LEARNER.pdf', width = 8, height = 8)
temp_target_learner <- svd_learner$v %*% t(svd_learner$v)
levelplot(myround(temp_target_learner, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Phenotype Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

pdf('PV-target-DLEARNER.pdf', width = 8, height = 8)
temp_target_dlearner <- svd_dlearner$v %*% t(svd_dlearner$v)
levelplot(myround(temp_target_dlearner, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Phenotype Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('D-LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()


################################################################################
## Plotting projection matrices onto the space of latent genotypic factors
################################################################################

min_val <- -0.00075; max_val <- 0.00075

pdf('PU-source.pdf', width = 8, height = 8)
levelplot(myround(PU_1_small, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('Source Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

pdf('PU-target.pdf', width = 8, height = 8)
levelplot(myround(PU_0_small, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('Target Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

pdf('PU-target-LEARNER.pdf', width = 8, height = 8)
levelplot(myround(PU_0_learner_small, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

pdf('PU-target-DLEARNER.pdf', width = 8, height = 8)
levelplot(myround(PU_0_dlearner_small, min = min_val, max = max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('D-LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()


################################################################################
## Plotting phenotype contribution scores
################################################################################

max_val <- 0.2
coul_1 <- rev(colorRampPalette(c("firebrick3", "white"))(20))
          
pdf('V2-source.pdf', width = 8, height = 8)
levelplot(pmin(svd1_chosen$v^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),  
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('Source Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('V2-target.pdf', width = 8, height = 8)
levelplot(pmin(svd0_chosen$v^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('Target Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('V2-target-LEARNER.pdf', width = 8, height = 8)
levelplot(pmin(svd_learner$v^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('V2-target-DLEARNER.pdf', width = 8, height = 8)
levelplot(pmin(svd_dlearner$v^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('D-LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()


################################################################################
## Plotting variant contribution scores
################################################################################

max_val <- 0.0005

pdf('U2-source.pdf', width = 8, height = 8)
levelplot(pmin(svd1_chosen$u^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('Source Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('U2-target.pdf', width = 8, height = 8)
levelplot(pmin(svd0_chosen$u^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('Target Population',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('U2-target-LEARNER.pdf', width = 8, height = 8)
levelplot(pmin(svd_learner$u^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()

pdf('U2-target-DLEARNER.pdf', width = 8, height = 8)
levelplot(pmin(svd_dlearner$u^2, max_val), 
          aspect="fill", 
          xlab = list(label = 'Variant Index', cex = cex.label),
          ylab = list(label = 'Component', cex = cex.label, at = 1:r_chosen),
          col.regions = coul_1, at = do.breaks(c(0,max_val),20), 
          main=list('D-LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales, at = 1:r_chosen)))
dev.off()


################################################################################
## Listing top phenotypes of the latent components
################################################################################

mycol <- colnames(theta_hat_0)
table_row <- function(dat, ind){
  abb <- mycol[order(dat$v[, ind]^2, decreasing = T)[1:5]]
  score <- round(sort(dat$v[, ind]^2, decreasing = T)[1:5], 2)
  res <- ""
  for (i in 1:length(abb)){
    res <- paste0(res, key[key$abb == abb[i], ]$Disease.name, 
                    ' (', 
                    score[i], 
                    '), ')
  }
  return(res)
}

table_row(svd0_chosen, 1)
table_row(svd0_chosen, 2)
table_row(svd0_chosen, 3)
table_row(svd0_chosen, 4)
table_row(svd0_chosen, 5)
table_row(svd0_chosen, 6)

table_row(svd1_chosen, 1)
table_row(svd1_chosen, 2)
table_row(svd1_chosen, 3)
table_row(svd1_chosen, 4)
table_row(svd1_chosen, 5)
table_row(svd1_chosen, 6)

table_row(svd_learner, 1)
table_row(svd_learner, 2)
table_row(svd_learner, 3)
table_row(svd_learner, 4)
table_row(svd_learner, 5)
table_row(svd_learner, 6)


table_row(svd_dlearner, 1)
table_row(svd_dlearner, 2)
table_row(svd_dlearner, 3)
table_row(svd_dlearner, 4)
table_row(svd_dlearner, 5)
table_row(svd_dlearner, 6)


################################################################################
## Listing top phenotypes of the latent components based on varimax rotation
################################################################################

mycol <- colnames(theta_hat_0)
table_row_loadings <- function(loadings, ind){
  abb <- mycol[order(loadings[, ind]^2, decreasing = T)[1:5]]
  score <- round(sort(loadings[, ind]^2, decreasing = T)[1:5], 2)
  res <- ""
  for (i in 1:length(abb)){
    res <- paste0(res, key[key$abb == abb[i], ]$Disease.name, 
                  ' (', 
                  score[i], 
                  '), ')
  }
  return(res)
}

rot_learner <- varimax(svd_learner$v)
rotated_loadings_learner <- rot_learner$loadings 
table_row_loadings(rotated_loadings_learner, 1)
table_row_loadings(rotated_loadings_learner, 2)
table_row_loadings(rotated_loadings_learner, 3)
table_row_loadings(rotated_loadings_learner, 4)
table_row_loadings(rotated_loadings_learner, 5)
table_row_loadings(rotated_loadings_learner, 6)

rot_dlearner <- varimax(svd_dlearner$v)
rotated_loadings_dlearner <- rot_dlearner$loadings 
table_row_loadings(rotated_loadings_dlearner, 1)
table_row_loadings(rotated_loadings_dlearner, 2)
table_row_loadings(rotated_loadings_dlearner, 3)
table_row_loadings(rotated_loadings_dlearner, 4)
table_row_loadings(rotated_loadings_dlearner, 5)
table_row_loadings(rotated_loadings_dlearner, 6)


################################################################################
## Listing top variants of the latent components
################################################################################

myrow <- rownames(theta_hat_0)
table_row_variant <- function(dat, ind){
  abb <- myrow[order(dat$u[, ind]^2, decreasing = T)[1:5]]
  score <- round(sort(dat$u[, ind]^2, decreasing = T)[1:5], 2)
  res <- ""
  for (i in 1:length(abb)){
    res <- paste0(res, abb[i], ' (', score[i], '), ')
  }
  return(res)
}
table_row_variant(svd0_chosen, 1)
table_row_variant(svd1_chosen, 1)
table_row_variant(svd_learner, 1)
table_row_variant(svd_dlearner, 1)

