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

plot_panel_unequal_rank <- function(data_tsvd, data_r4, data_r8, title, ylim) {
  plot(var_ratio, data_tsvd$tsvd, col = cols[1], ylim = ylim,
       xlab = expression(sigma[0]^2/sigma[1]^2),
       ylab = 'Estimation error', main = title,
       type = 'b', pch = 16, lwd = 2, lty = 1,
       cex.axis = 1.18, cex.lab = 1.18, cex.main = 1.25)
  grid(col = "gray85", lty = 1)
  points(var_ratio, data_r4$dlearner, col = cols[2], type = 'b', pch = 16, lwd = 2, lty = 2)
  points(var_ratio, data_r8$dlearner, col = cols[2], type = 'b', pch = 2, lwd = 2, lty = 2)
  points(var_ratio, data_r4$learner, col = cols[3], type = 'b', pch = 16, lwd = 2, lty = 6)
  points(var_ratio, data_r8$learner, col = cols[3], type = 'b', pch = 2, lwd = 2, lty = 6)
}


################################################################################
## Independent noise simulations
################################################################################

load('simres-same.RData')
load('simres-dif-moderate.RData')
load('simres-dif-verylow.RData')
var_ratio <- var_0_all / var_1_all

pdf(paste0('Sim-Res-Rectangle.pdf'), height = 7.65, width = 5.1)
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
## Independent noise simulations: Square matrix
################################################################################

load('simres-same-square.RData')
load('simres-dif-moderate-square.RData')
load('simres-dif-verylow-square.RData')
var_ratio <- var_0_all / var_1_all

pdf(paste0('Sim-Res-Square.pdf'), height = 7.65, width = 5.1)
ylim <- c(0, 55)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow = 4, byrow = TRUE), heights = c(1, 1, 1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.5, 1, 0))

panel_plot(res_1_same_square, "High Similarity (Rank = 4)", ylim)
panel_plot(res_1_same_square, "High Similarity (Rank = 8)", ylim)
panel_plot(res_1_dif_moderate_square, "Moderate Similarity (Rank = 4)", ylim)
panel_plot(res_2_dif_moderate_square, "Moderate Similarity (Rank = 8)", ylim)
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


################################################################################
## Unequal-rank simulations
################################################################################

load('simres-unequal-rank.RData')
var_ratio <- var_0_all / var_1_all
ylim <- c(0, max(c(res_4_8$r4$tsvd, res_4_8$r4$dlearner, res_4_8$r4$learner,
                   res_4_8$r8$tsvd, res_4_8$r8$dlearner, res_4_8$r8$learner,
                   res_8_4$r4$tsvd, res_8_4$r4$dlearner, res_8_4$r4$learner,
                   res_8_4$r8$tsvd, res_8_4$r8$dlearner, res_8_4$r8$learner), 
                 na.rm = TRUE) * 1.05)

pdf(paste0('Sim-Res-Unequal-Rank.pdf'), height = 4.25, width = 7)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.5, 1, 0))

plot_panel_unequal_rank(res_4_8$r4, res_4_8$r4, res_4_8$r8, "Target Rank = 4\nSource Rank = 8", ylim)
plot_panel_unequal_rank(res_8_4$r8, res_8_4$r4, res_8_4$r8, "Target Rank = 8\nSource Rank = 4", ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Target-Only SVD", "D-LEARNER (rank = 4)", "D-LEARNER (rank = 8)",
                            "LEARNER (rank = 4)", "LEARNER (rank = 8)"),
       col = c(cols[1], cols[2], cols[2], cols[3], cols[3]), lty = c(1, 2, 2, 6, 6),
       pch = c(16, 16, 2, 16, 2), bty = "n", ncol = 2, cex = 1.075)
dev.off()


################################################################################
## Multi-source simulations (independent noise, rank 4, high similarity)
################################################################################
library(RColorBrewer)
cols <- brewer.pal(3, "Set1")[c(1,3,2)]

load('simres-multisource.RData')
res_multisource <- res_indep_r4_high_0.1

# --- data ---
err <- c(
  res_multisource$tsvd,
  res_multisource$dlearner_a,    res_multisource$dlearner_b,
  res_multisource$dlearner_list, res_multisource$dlearner_list_w,
  res_multisource$learner_a,     res_multisource$learner_b,
  res_multisource$learner_list,  res_multisource$learner_list_w
)

bar_cols <- c(cols[1], rep(cols[2], 4), rep(cols[3], 4))
pchs     <- c(18, 16, 2, 15, 17, 16, 2, 15, 17)
n        <- length(err)
xs       <- seq_len(n)
ylim     <- c(0, max(err, na.rm = TRUE) * 1.12)

xlabs <- c("Target-Only SVD",
           "Source A", "Source B", "Multisource (unweighted)", "Multisource (weighted)",
           "Source A", "Source B", "Multisource (unweighted)", "Multisource (weighted)")

pdf('Sim-Res-Multisource.pdf', height = 4.5, width = 7)
layout(matrix(c(1), nrow = 1))
par(mar = c(9, 4, 1, 2), mgp = c(2.5, 1, 0))

# --- main panel ---
plot(NA, xlim = c(0.5, n + 0.5), ylim = ylim,
     xaxt = 'n', xlab = '', ylab = 'Estimation error',
     cex.axis = 1.07, cex.lab = 1.1, cex.main = 1.1)

grid(nx = NA, ny = NULL, col = 'gray85', lty = 1)

# subtle group shading
rect(0.5, 0, 1.5, ylim[2] * 2, col = adjustcolor(cols[1], alpha.f = 0.075), border = NA)
rect(1.5, 0, 5.5, ylim[2] * 2, col = adjustcolor(cols[2], alpha.f = 0.075), border = NA)
rect(5.5, 0, 9.5, ylim[2] * 2, col = adjustcolor(cols[3], alpha.f = 0.075), border = NA)

# vertical dividers
abline(v = c(1.5, 5.5), col = 'gray70', lty = 3, lwd = 1.2)

# connecting lines within groups
lines(2:5, err[2:5], col = cols[2], lwd = 1.5, lty = 2)
lines(6:9, err[6:9], col = cols[3], lwd = 1.5, lty = 2)

# points
points(xs, err, col = bar_cols, pch = pchs, cex = 1.6, lwd = 1.5)

# x-axis with rotated labels
axis(1, at = xs, labels = FALSE)
text(xs, par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
     labels = xlabs, srt = 35, adj = 1, xpd = TRUE, cex = 0.8)

# group labels below x-axis
mtext("SVD",       side = 1, at = 1,          line = 6.5,
      cex = 1.05, col = cols[1], font = 2)
mtext("D-LEARNER", side = 1, at = mean(2:5), line = 6.5,
      cex = 1.05, col = cols[2], font = 2)
mtext("LEARNER",   side = 1, at = mean(6:9), line = 6.5,
      cex = 1.05, col = cols[3], font = 2)


dev.off()


################################################################################
## Multi-source high vs low similarity (5000 x 50): A same latent space, B dif-verylow scale
################################################################################
library(RColorBrewer)
cols <- brewer.pal(3, "Set1")[c(1, 3, 2)]

load('simres-multisource-high-low.RData')
res_hl <- res_multisource_high_low

err_hl <- c(
  res_hl$tsvd,
  res_hl$dlearner_a,    res_hl$dlearner_b,
  res_hl$dlearner_list, res_hl$dlearner_list_w,
  res_hl$learner_a,     res_hl$learner_b,
  res_hl$learner_list,  res_hl$learner_list_w
)

bar_cols <- c(cols[1], rep(cols[2], 4), rep(cols[3], 4))
pchs     <- c(18, 16, 2, 15, 17, 16, 2, 15, 17)
n_hl     <- length(err_hl)
xs_hl    <- seq_len(n_hl)
ylim_hl  <- c(0, max(err_hl, na.rm = TRUE) * 1.12)

xlabs_hl <- c("Target-Only SVD",
              "Source A", "Source B", "Multisource (unweighted)", "Multisource (weighted)",
              "Source A", "Source B", "Multisource (unweighted)", "Multisource (weighted)")

pdf('Sim-Res-Multisource-High-Low.pdf', height = 4.5, width = 7)
layout(matrix(c(1), nrow = 1))
par(mar = c(9, 4, 2.2, 2), mgp = c(2.5, 1, 0))

plot(NA, xlim = c(0.5, n_hl + 0.5), ylim = ylim_hl,
     xaxt = 'n', xlab = '', ylab = 'Estimation error',
     cex.axis = 1.07, cex.lab = 1.1, cex.main = 1.1)

grid(nx = NA, ny = NULL, col = 'gray85', lty = 1)

rect(0.5, 0, 1.5, ylim_hl[2] * 2, col = adjustcolor(cols[1], alpha.f = 0.075), border = NA)
rect(1.5, 0, 5.5, ylim_hl[2] * 2, col = adjustcolor(cols[2], alpha.f = 0.075), border = NA)
rect(5.5, 0, 9.5, ylim_hl[2] * 2, col = adjustcolor(cols[3], alpha.f = 0.075), border = NA)

abline(v = c(1.5, 5.5), col = 'gray70', lty = 3, lwd = 1.2)

lines(2:5, err_hl[2:5], col = cols[2], lwd = 1.5, lty = 2)
lines(6:9, err_hl[6:9], col = cols[3], lwd = 1.5, lty = 2)

points(xs_hl, err_hl, col = bar_cols, pch = pchs, cex = 1.6, lwd = 1.5)

axis(1, at = xs_hl, labels = FALSE)
text(xs_hl, par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
     labels = xlabs_hl, srt = 35, adj = 1, xpd = TRUE, cex = 0.8)

mtext("SVD",       side = 1, at = 1,          line = 6.5,
      cex = 1.05, col = cols[1], font = 2)
mtext("D-LEARNER", side = 1, at = mean(2:5), line = 6.5,
      cex = 1.05, col = cols[2], font = 2)
mtext("LEARNER",   side = 1, at = mean(6:9), line = 6.5,
      cex = 1.05, col = cols[3], font = 2)

dev.off()


################################################################################
## Independent noise: U high similarity / V perturbed (row/col penalties, learnerv2)
################################################################################
library(RColorBrewer)
cols <- brewer.pal(4, "Set1")[c(1, 3, 2, 4)]

plot_panel_uhigh_flex <- function(data, title, ylim) {
  plot(var_ratio, data$tsvd, col = cols[1], ylim = ylim,
       xlab = expression(sigma[0]^2/sigma[1]^2),
       ylab = 'Estimation error', main = title,
       type = 'b', pch = 16, lwd = 2, lty = 1,
       cex.axis = 1.18, cex.lab = 1.18, cex.main = 1.25)
  grid(col = "gray85", lty = 1)
  points(var_ratio, data$dlearner, col = cols[2], type = 'b', pch = 16, lwd = 2, lty = 2)
  points(var_ratio, data$learner_tied, col = cols[3], type = 'b', pch = 16, lwd = 2, lty = 6)
  points(var_ratio, data$learner_sep, col = cols[4], type = 'b', pch = 16, lwd = 2, lty = 6)
}

load('simres-uhigh-vlow.RData')
load('simres-uhigh-vmoderate.RData')
var_ratio <- res_uhigh_vlow$var_0_all / res_uhigh_vlow$var_1_all

ylim <- c(0, max(c(res_uhigh_vlow$tsvd, res_uhigh_vlow$dlearner,
                   res_uhigh_vlow$learner_tied, res_uhigh_vlow$learner_sep,
                   res_uhigh_vmoderate$tsvd, res_uhigh_vmoderate$dlearner,
                   res_uhigh_vmoderate$learner_tied, res_uhigh_vmoderate$learner_sep),
                 na.rm = TRUE) * 1.05)

pdf('Sim-Res-Uhigh-V.pdf', height = 4.25, width = 7)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.5, 1, 0))

plot_panel_uhigh_flex(res_uhigh_vmoderate, "High Similarity in U\nModerate Similarity in V", ylim)
plot_panel_uhigh_flex(res_uhigh_vlow, "High Similarity in U\nLow Similarity in V", ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("Target-Only SVD", "D-LEARNER", "LEARNER (same penalty)", "LEARNER (different penalty)"),
       col = cols, lty = c(1, 2, 6, 6), pch = c(16, 16, 16, 16),
       bty = "n", ncol = 2, cex = 1.1)
dev.off()


################################################################################
## Rectangular (5000×50): high vs moderate latent similarity (learnerv2 tied / separate)
################################################################################
library(RColorBrewer)
cols <- brewer.pal(4, "Set1")[c(1, 3, 2, 4)]

plot_panel_rect_flex <- function(data, title, ylim) {
  plot(var_ratio, data$tsvd, col = cols[1], ylim = ylim,
       xlab = expression(sigma[0]^2/sigma[1]^2),
       ylab = 'Estimation error', main = title,
       type = 'b', pch = 16, lwd = 2, lty = 1,
       cex.axis = 1.18, cex.lab = 1.18, cex.main = 1.25)
  grid(col = "gray85", lty = 1)
  points(var_ratio, data$dlearner, col = cols[2], type = 'b', pch = 16, lwd = 2, lty = 2)
  points(var_ratio, data$learner_tied, col = cols[3], type = 'b', pch = 16, lwd = 2, lty = 6)
  points(var_ratio, data$learner_sep, col = cols[4], type = 'b', pch = 16, lwd = 2, lty = 6)
}

load('simres-rect-flex-high.RData')
load('simres-rect-flex-moderate.RData')
var_ratio <- res_rect_high$var_0_all / res_rect_high$var_1_all

ylim <- c(0, max(c(res_rect_high$tsvd, res_rect_high$dlearner,
                   res_rect_high$learner_tied, res_rect_high$learner_sep,
                   res_rect_moderate$tsvd, res_rect_moderate$dlearner,
                   res_rect_moderate$learner_tied, res_rect_moderate$learner_sep),
                 na.rm = TRUE) * 1.05)

pdf('Sim-Res-Rect-Flex.pdf', height = 4.25, width = 7)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))
par(mar = c(4, 4, 3.5, 2), mgp = c(2.5, 1, 0))

plot_panel_rect_flex(res_rect_high, "High Similarity", ylim)
plot_panel_rect_flex(res_rect_moderate, "Moderate Similarity", ylim)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("Target-Only SVD", "D-LEARNER", "LEARNER (same penalty)", "LEARNER (different penalty)"),
       col = cols, lty = c(1, 2, 6, 6), pch = c(16, 16, 16, 16),
       bty = "n", ncol = 2, cex = 1.1)
dev.off()