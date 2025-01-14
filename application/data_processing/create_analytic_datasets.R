rm(list = ls())

theta_hat_0 <- read.delim('filtered_table_BBJ.txt', header = TRUE)
theta_hat_1 <- read.delim('filtered_table_EUR.txt', header = TRUE)

theta_hat_0_app <- read.delim('filtered_table_append_BBJ.txt', header = TRUE)
theta_hat_1_app <- read.delim('filtered_table_append_EUR.txt', header = TRUE)

theta_hat_1 <- theta_hat_1[theta_hat_1$X %in% theta_hat_1_app$X,]
theta_hat_0 <- theta_hat_0[theta_hat_0$X %in% theta_hat_1_app$X,]
theta_hat_0_app <- theta_hat_0_app[theta_hat_0_app$X %in% theta_hat_1_app$X,]

theta_hat_0 <- cbind(theta_hat_0, theta_hat_0_app[, -1])
theta_hat_1 <- cbind(theta_hat_1, theta_hat_1_app[, -1])

snps_0 <- theta_hat_0[, 1]
snps_1 <- theta_hat_1[, 1]

theta_hat_0 <- theta_hat_0[, -1] 
theta_hat_1 <- theta_hat_1[, -1] 

rownames(theta_hat_0) <- snps_0
rownames(theta_hat_1) <- snps_1

chr_0 <- as.numeric(sub(":.*", "", snps_0))
pos_0 <- as.numeric(sub(".*:", "", snps_0))
order_0 <- order(chr_0, pos_0)

chr_1 <- as.numeric(sub(":.*", "", snps_1))
pos_1 <- as.numeric(sub(".*:", "", snps_1))
order_1 <- order(chr_1, pos_1)

theta_hat_0 <- theta_hat_0[, order(colnames(theta_hat_0))]
theta_hat_0 <- theta_hat_0[order_0, ]

theta_hat_1 <- theta_hat_1[, order(colnames(theta_hat_1))]
theta_hat_1 <- theta_hat_1[order_1, ]

theta_hat_0[is.na(theta_hat_0)] <- 0

# Newly added this
theta_hat_0 <- as.matrix(theta_hat_0)
theta_hat_1 <- as.matrix(theta_hat_1) 

p <- nrow(theta_hat_0)
q <- ncol(theta_hat_0)
r_chosen <- 8

set.seed(1234)
heldout_indices <- sample(1:(p * q), size = p * q, replace = FALSE)
heldout_set <- vector(mode = 'list', length = 5)
heldout_set[[1]] <- heldout_indices[1:floor(p * q / 5)]
heldout_set[[2]] <- heldout_indices[floor(p * q / 5 + 1):floor(2 * p * q / 5)]
heldout_set[[3]] <- heldout_indices[floor(2 * p * q / 5 + 1):floor(3 * p * q / 5)]
heldout_set[[4]] <- heldout_indices[floor(3 * p * q / 5 + 1):floor(4 * p * q / 5)]
heldout_set[[5]] <- heldout_indices[floor(4 * p * q / 5 + 1):(p * q)]

save.image('../dat.RData')
