## Plotting the lambda selection
rm(list = ls())
load('lambda_small.RData')
load('../dat.RData')
library('lattice')
library('RColorBrewer')

n_lambdas <- 10
log_lambdas1 <- round(seq(from = 1, to = 4, length.out = n_lambdas), 2)
log_lambdas2 <- round(seq(from = -2, to = 1, length.out = n_lambdas), 2)
min_ind <- which(out_small_lambda$res == min(out_small_lambda$res, na.rm = TRUE), arr.ind = TRUE)

coul <- rev(colorRampPalette(brewer.pal(9, "Reds"))(100))

pdf('../Lambdas-4fold.pdf', width = 7, height = 6)
levelplot(out_small_lambda$res / (p * q), aspect = 'fill',
          col.regions = coul, 
          xlab = expression(log(lambda[1])), ylab = expression(log(lambda[2])), 
          alpha.regions = 0.85, 
          scales = list(x = list(labels = log_lambdas1, at = 1:n_lambdas), 
                        y = list(labels = log_lambdas2, at = 1:n_lambdas)))

trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(min_ind[1], min_ind[2], pch=4, cex=2, col = 'white')
trellis.unfocus()
dev.off()

## Plotting convergence of the optimization algorithm
rm(list = ls())
load('learner.RData')

pdf('../ObjFunc.pdf', width = 6, height = 4)
par(mar = c(5.1, 4.6, 4.1, 2.1))
plot(fit$objective_values[1:100] / 100000, type = 'l', col = 'blue', lwd = 2, 
     xlab = 'Iteration number', ylab = expression(paste("Objective function (10"^6, ")")))
dev.off()

## Plotting LEARNER and D-LEARNER estimates
rm(list = ls())
load('../dat.RData')
load('learner.RData')
load('DLEARNER-estimate-BBJ.RData')
load('../key.RData')

fit$learner_estimate <- fit$learner_estimate[, order(key$ICD_cat)]
theta_hat_0_dlearner <- theta_hat_0_dlearner[, order(key$ICD_cat)]

coul <- c(colorRampPalette(c("firebrick3", "white"))(10),
          colorRampPalette(c("white", "dodgerblue3"))(10))
cex.main <- 2.25; cex.label <- 1.75; cex.scales <- 1.5

set.seed(1234)
myind <- sort(sample.int(p, 1000, replace = F))

temp <- t(fit$learner_estimate[myind,])
colnames(temp) <- rownames(temp) <- NULL
min_val <- -1.25; max_val <- 1.25
pdf('../Thetahat-LEARNER.pdf', width = 8, height = 8)
levelplot(temp, aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()

temp <- t(theta_hat_0_dlearner[myind,])
colnames(temp) <- rownames(temp) <- NULL
min_val <- -0.2; max_val <- 0.2
pdf('../Thetahat-DLEARNER.pdf', width = 8, height = 8)
levelplot(temp, aspect="fill", 
          xlab = list(label = 'Phenotype Index', cex = cex.label),
          ylab = list(label = 'Variant Index', cex = cex.label),  
          col.regions = coul, at = do.breaks(c(min_val,max_val),20), 
          main=list('D-LEARNER Estimate',cex = cex.main), 
          scales = list(x=list(cex=cex.scales),y=list(cex=cex.scales)))
dev.off()