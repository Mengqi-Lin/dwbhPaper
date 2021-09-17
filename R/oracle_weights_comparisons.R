Pi1 <- rbind(c(rep(0.1, 5), rep(0.01, 5)), c(rep(0.1, 5), rep(0.01, 5)))
MU1 <- rbind(rep(1, 10), c(seq(2,6,1), seq(2,6,1)))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 10)
gw <- matrix(NA, 2, 10)
for(i in 1:10) {
  pi0_g = 1-Pi1[,i]
  pi1 = sum(Pi1[,i])*0.5
  mu1_g = MU1[,i]
  ow[,i] <- unique(oracle.weights(alpha, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g))
  gw[,i] <- (1-pi0_g)/pi0_g/pi1
} 

FDR.filename <- "../figs/oracle_weights1.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
plot(c(0, 3), c(0, 2.5), type = 'n',
     xlab = "group1",
     ylab = "group2",
     main = expression(pi[1]~"fixed, " ~ mu[1]~"varied, "~ alpha~ "= 0.05"), axes = T)
for(i in 1:10) {
  points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
  text(ow[,i][1] + 0.15, ow[,i][2], labels = i)
}
abline(coef = c(0,1))
dev.off()










Pi1 <- rbind(seq(0.01, 0.2, length.out = 10), seq(0.01, 0.2, length.out = 10))
MU1 <- rbind(rep(2, 10), rep(4, 10))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 10)
gw <- matrix(NA, 2, 10)
for(i in 1:10) {
  pi0_g = 1-Pi1[,i]
  pi1 = sum(Pi1[,i])*0.5
  mu1_g = MU1[,i]
  ow[,i] <- unique(oracle.weights(alpha, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g))
  gw[,i] <- (1-pi0_g)/pi0_g/pi1
} 

plot(c(0, 3), c(0, 2.5), type = 'n',
     xlab = "group1",
     ylab = "group2",
     main = expression(alpha~ "= 0.05"), axes = T)
for(i in 1:10) {
  points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
  text(ow[,i][1] + 0.15, ow[,i][2], labels = i)
}
abline(coef = c(0,1))