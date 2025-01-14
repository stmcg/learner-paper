rm(list = ls())

library('ScreeNOT')

load('dat.RData')
res_0 <- adaptiveHardThresholding(Y = theta_hat_0, k = 30) # rank 6
res_1 <- adaptiveHardThresholding(Y = theta_hat_1, k = 30) # rank 8

svd_0 <- svd(theta_hat_0)
svd_1 <- svd(theta_hat_1)

plot(1:length(svd_0$d), svd_0$d, 
     xlab = 'Component number', ylab = 'Singular value')
abline(h = res_0$Topt)

pdf('ScreePlot.pdf', width = 6, height = 5)
plot(1:length(svd_1$d), svd_1$d, 
     xlab = 'Component number', ylab = 'Singular value',
     pch = 19, cex = 0.5, type = 'o')
abline(h = res_1$Topt, lty = 2, lwd = 2, col = 'red')
text(120, res_1$Topt+20, 'ScreeNOT threshold', col = 'red')
dev.off()