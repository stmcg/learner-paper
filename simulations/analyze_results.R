rm(list = ls())
library(RColorBrewer)
cols <- brewer.pal(3, "Set1")[c(1,3,2)]

panel_plot <- function(data, title, ylim) {
  plot(var_ratio, data$tsvd, col = cols[1], ylim = ylim, 
       xlab = expression(sigma[0]^2/sigma[1]^2), 
       ylab = 'Estimation error', main = title, 
       type = 'b', pch = 16, lwd = 2, lty = 1, 
       cex.axis = 1.07, cex.lab = 1.1, cex.main = 1.1)
  grid(col = "gray85", lty = 1)
  points(var_ratio, data$dlearner, col = cols[2], type = 'b', pch = 16, lwd = 2, lty = 2)
  points(var_ratio, data$learner, col = cols[3], type = 'b', pch = 16, lwd = 2, lty = 6)
}

plot_panel_cor <- function(data, title, ylim) {
  plot(var_ratio, data$tsvd, col = cols[1], ylim = ylim, 
       xlab = expression(sigma[0]^2/sigma[1]^2), 
       ylab = 'Estimation error', main = title, 
       type = 'b', pch = 16, lwd = 2, lty = 1, 
       cex.axis = 1.325, cex.lab = 1.325, cex.main = 1.45)
  grid(col = "gray85", lty = 1)
  segments(var_ratio, data$tsvd - data$tsvd_sd, var_ratio, data$tsvd + data$tsvd_sd, col = cols[1])
  points(var_ratio, data$dlearner, col = cols[2], type = 'b', pch = 16, lwd = 2, lty = 2)
  segments(var_ratio, data$dlearner - data$dlearner_sd, var_ratio, data$dlearner + data$dlearner_sd, col = cols[2])
  points(var_ratio, data$learner, col = cols[3], type = 'b', pch = 16, lwd = 2, lty = 6)
  segments(var_ratio, data$learner - data$learner_sd, var_ratio, data$learner + data$learner_sd, col = cols[3])
}


################################################################################
## Independent noise simulations
################################################################################

load('simres-same.RData')
load('simres-dif-moderate.RData')
load('simres-dif-verylow.RData')
var_ratio <- var_0_all / var_1_all

pdf(paste0('Sim-Res-Independent.pdf'), height = 7.65, width = 5.1)
ylim <- c(0, 100)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow = 4, byrow = TRUE), heights = c(1, 1, 1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.5, 1, 0))

panel_plot(res_1_same, "High Similarity (Rank = 4)", ylim)
panel_plot(res_2_same, "High Similarity (Rank = 8)", ylim)
panel_plot(res_1_dif_moderate, "Moderate Similarity (Rank = 4)", ylim)
panel_plot(res_2_dif_moderate, "Moderate Similarity (Rank = 8)", ylim)
panel_plot(res_1_dif_verylow, "Low Similarity (Rank = 4)", ylim)
panel_plot(res_2_dif_verylow, "Low Similarity (Rank = 8)", ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Target-Only SVD", "D-LEARNER", "LEARNER"), col = cols, lty = c(1, 2, 6), 
       pch = 16, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()



################################################################################
## Correlated noise simulations
################################################################################

load('simres-same-cor.RData')
load('simres-same-cor2.RData')
load('simres-same-cor3.RData')

load('simres-dif-moderate-cor.RData')
load('simres-dif-moderate-cor2.RData')
load('simres-dif-moderate-cor3.RData')

load('simres-dif-verylow-cor.RData')
load('simres-dif-verylow-cor2.RData')
load('simres-dif-verylow-cor3.RData')


pdf(paste0('Sim-Res-Correlation.pdf'), height = 8.5, width = 8)
var_ratio <- var_0_all / var_1_all
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10), nrow = 4, byrow = TRUE), heights = c(1, 1, 1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.8, 1, 0))

ylim <- c(0, 150)
plot_panel_cor(res_same_cor1, expression(bold("High Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_same_cor2, expression(bold("High Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_same_cor3, expression(bold("High Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

plot_panel_cor(res_dif_moderate_cor1, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_dif_moderate_cor2, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_dif_moderate_cor3, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

plot_panel_cor(res_dif_verylow_cor1, expression(bold("Low Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_dif_verylow_cor2, expression(bold("Low Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_dif_verylow_cor3, expression(bold("Low Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Target-Only SVD", "D-LEARNER", "LEARNER"), col = cols, lty = c(1, 2, 6), 
       pch = 16, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()



################################################################################
## Correlated noise simulations - external data set
################################################################################

load('simres-same-cor-ext.RData')
load('simres-same-cor2-ext.RData')
load('simres-same-cor3-ext.RData')

load('simres-dif-moderate-cor-ext.RData')
load('simres-dif-moderate-cor2-ext.RData')
load('simres-dif-moderate-cor3-ext.RData')

load('simres-dif-verylow-cor-ext.RData')
load('simres-dif-verylow-cor2-ext.RData')
load('simres-dif-verylow-cor3-ext.RData')


pdf(paste0('Sim-Res-Correlation-External.pdf'), height = 8.5, width = 8)
var_ratio <- var_0_all / var_1_all
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10), nrow = 4, byrow = TRUE), heights = c(1, 1, 1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.8, 1, 0))

ylim <- c(0, 150)
plot_panel_cor(res_same_cor1_ext, expression(bold("High Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_same_cor2_ext, expression(bold("High Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_same_cor3_ext, expression(bold("High Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

plot_panel_cor(res_dif_moderate_cor1_ext, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_dif_moderate_cor2_ext, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_dif_moderate_cor3_ext, expression(bold("Moderate Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

plot_panel_cor(res_dif_verylow_cor1_ext, expression(bold("Low Similarity")~"("*bold(rho)*"=0.1)"), ylim = ylim)
plot_panel_cor(res_dif_verylow_cor2_ext, expression(bold("Low Similarity")~"("*bold(rho)*"=0.25)"), ylim = ylim)
plot_panel_cor(res_dif_verylow_cor3_ext, expression(bold("Low Similarity")~"("*bold(rho)*"=0.5)"), ylim = ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Target-Only SVD", "D-LEARNER", "LEARNER"), col = cols, lty = c(1, 2, 6), 
       pch = 16, bty = "n", horiz = TRUE, cex = 1.2)
dev.off()
