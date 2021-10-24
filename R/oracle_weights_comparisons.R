## Combination 1: Same pi1 in two groups, vary on mu1.

Pi1 <- rbind(c(rep(0.1, 3), rep(0.05, 3)), c(rep(0.1, 3), rep(0.05, 3)))
MU1 <- rbind(rep(2, 6), c(seq(3,5,1), seq(3,5,1)))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 6)
gw <- matrix(NA, 2, 6)


FDR.filename <- "../figs/CombinationWeights_same_pi1.pdf"
pdf(FDR.filename, width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
for(alpha in c(0.05, 0.1, 0.2)) {
  for(i in 1:6) {
    pi0_g = 1-Pi1[,i]
    pi1 = sum(Pi1[,i])*0.5
    mu1_g = MU1[,i]
    ow[,i] <- unique(oracle.weights(alpha, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g, F))
    gw[,i] <- (1-pi0_g)/pi0_g/pi1
  } 
  plot(c(0, 2), c(0, 2), type = 'n',
       xlab = "group1",
       ylab = "group2",
       main = bquote(pi[1]~"fixed, " ~ mu[1]~"varied, "~ alpha~" = " ~.(alpha)), axes = T)
  for(i in 1:6) {
    points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
    text(ow[,i][1] + 0.15, ow[,i][2], labels = i)
  }
  abline(coef = c(0,1))
  legend("topright", legend = paste(Pi1[1,],  " vs ", Pi1[2], ", ", MU1[1,], " vs ",MU1[2,]))
}
dev.off()

comb <- matrix
dist = c()
for(i in 1:6) {
  v = gw[, i] - ow[,i]
  dist[i] = v%*%v
}
ind <- which(dist > 1)
temp <- rbind(Pi1[,ind], MU1[,ind])
comb <- temp

## Combination 2: Same mu1 in two groups, vary on pi1.

Pi1 <- rbind(rep(0.05, 10), rep(seq(0.1, 0.2, length.out = 5), 2))
MU1 <- rbind(c(rep(2,5), rep(3,5)), c(rep(2,5), rep(3,5)))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 10)
gw <- matrix(NA, 2, 10)

FDR.filename <- "../figs/CombinationWeights_same_mu1.pdf"
pdf(FDR.filename, width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
for(alpha in c(0.05, 0.1, 0.2)) {
  for(i in 1:10) {
    pi0_g = 1-Pi1[,i]
    pi1 = sum(Pi1[,i])*0.5
    mu1_g = MU1[,i]
    ow[,i] <- unique(oracle.weights(alpha = 0.05, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g, F))
    gw[,i] <- (1-pi0_g)/pi0_g/pi1
  } 
  plot(c(0, 3), c(0, 2.5), type = 'n',
       xlab = "group1",
       ylab = "group2",
       main = bquote(mu[1]~ "fixed, " ~ pi[1]~"varied, "~ alpha~" = " ~.(alpha)), axes = T)
  for(i in 1:10) {
    points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
    text(ow[,i][1] + 0.04, ow[,i][2], labels = i)
  }
  legend("topright", legend = paste(MU1[1,], ", ", Pi1[1,], " vs ",Pi1[2,]))
  abline(coef = c(0,1))
}
dev.off()


dist = c()
for(i in 1:10) {
  v = gw[, i] - ow[,i]
  dist[i] = v%*%v
}
ind <- which(dist > 1)
temp <- rbind(Pi1[,ind], MU1[,ind])
comb <- cbind(comb, temp)

## Combination 3: mu1 increasing in pi1.
Pi1 <- rbind(c(rep(0.01, 4), rep(0.05, 4)), c(seq(0.05, 0.2, length.out = 4), seq(0.05, 0.2, length.out = 4)))
MU1 <- rbind(rep(1, 8), c(seq(2,5,1), seq(2,5,1)))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 8)
gw <- matrix(NA, 2, 8)
FDR.filename <- "../figs/CombinationWeights_inc.pdf"
pdf(FDR.filename, width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
for(alpha in c(0.05, 0.1, 0.25)) {
  for(i in 1:8) {
    pi0_g = 1-Pi1[,i]
    pi1 = sum(Pi1[,i])*0.5
    mu1_g = MU1[,i]
    ow[,i] <- unique(oracle.weights(alpha, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g, F))
    gw[,i] <- (1-pi0_g)/pi0_g/pi1*(1-pi1)
  } 
  plot(c(0, 2.5), c(0, 2.5), type = 'n',
       xlab = "group1",
       ylab = "group2",
       main = bquote(mu[1]~ " decreasing in " ~ pi[1]~" "~ alpha~" = " ~.(alpha)), axes = T)
  for(i in 1:8) {
    points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
    text(ow[,i][1] + 0.1, ow[,i][2], labels = i)
  }
  legend("topright", legend = paste0(MU1[1,], " vs ", MU1[2,], 
                                     ", ",
                                     Pi1[1,], " vs ",Pi1[2,]))
  abline(coef = c(0,1))
}

dev.off()


dist = c()
for(i in 1:8) {
  v = gw[, i] - ow[,i]
  dist[i] = v%*%v
}
ind <- which(dist > 1)
temp <- rbind(Pi1[,ind], MU1[,ind])
comb <- cbind(comb, temp)

## Combination 4: mu1 decreasing in pi1.

Pi1 <- rbind(c(seq(0.05, 0.2, length.out = 4), seq(0.05, 0.2, length.out = 4)), c(rep(0.01, 4), rep(0.05, 4)))
MU1 <- rbind(rep(1, 8), c(rev(seq(2,5,1)), rev(seq(2,5,1))))
n_g <- c(1000, 1000)
ow <- matrix(NA, 2, 8)
gw <- matrix(NA, 2, 8)
FDR.filename <- "../figs/CombinationWeights_dec.pdf"
pdf(FDR.filename, width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
for(alpha in c(0.05, 0.1, 0.25)) {
  for(i in 1:8) {
    pi0_g = 1-Pi1[,i]
    pi1 = sum(Pi1[,i])*0.5
    mu1_g = MU1[,i]
    ow[,i] <- unique(oracle.weights(alpha, n_g = n_g, pi0_g = pi0_g, mu1_g = mu1_g, F))
    gw[,i] <- (1-pi0_g)/pi0_g/pi1*(1-pi1)
  } 
  plot(c(0, 2.5), c(0, 2.5), type = 'n',
       xlab = "group1",
       ylab = "group2",
       main = bquote(mu[1]~ " decreasing in " ~ pi[1]~" "~ alpha~" = " ~.(alpha)), axes = T)
  for(i in 1:8) {
    points(c(ow[,i][1], gw[,i][1]), c(ow[,i][2], gw[,i][2]), col = c("red", "black"), type = "b")
    text(ow[,i][1] + 0.1, ow[,i][2], labels = i)
  }
  legend("topright", legend = paste0(MU1[1,], " vs ", MU1[2,], 
                                     ", ",
                                     Pi1[1,], " vs ",Pi1[2,]))
  abline(coef = c(0,1))
}

dev.off()



dist = c()
for(i in 1:8) {
  v = gw[, i] - ow[,i]
  dist[i] = v%*%v
}
ind <- which(dist > 1)
if(length(ind)==1){
  temp <- c(Pi1[,ind], MU1[,ind])
}

comb <- cbind(comb, temp)
