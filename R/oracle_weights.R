
oracle.weights_mvgauss <- function(alpha, n_g, pi0_g, mu1_g, pi0Est = T, side) {
  pi_g <- n_g/sum(n_g)
  #lambda_alpha <- solve.mFDR(alpha, n_g, pi0_g, mu1_g)
  mFDR <- function(c) {
    V <- c()
    R <- c()
    for(i in 1:length(pi0_g)) {
      pi0 = pi0_g[i]
      mu1 = mu1_g[i]
      n = n_g[i]
      t = lfdr.inverse_mvgauss(c, pi0, mu1, side = side)
      V[i] = n*pi0*t
      R[i] = V[i] + n*(1-pi0)*nonnull.cdf_mvgauss(t, mu1)
    }
    return(sum(V)/sum(R))
  }
  lambda_alpha <- vecbinsolv(zf = alpha, fun = mFDR , tlo = 0, thi = 1, nits = 30)
  w <- c()
  t <- c()
  for(i in 1:length(n_g)) {
    t[i] <- lfdr.inverse_mvgauss(lambda_alpha, pi0 = pi0_g[i], mu1 = mu1_g[i], side = side)
  }
  if(pi0Est) {
    t <- t/sum(t*pi0_g*pi_g)
  } else{
    t <- t/mean(t)
  }
  
  for(i in 1:length(n_g)) {
    w <- c(w, rep(t[i], n_g[i]))
  }
  return(lambda_alpha, )
}

lfdr.inverse_mvgauss <- function(lfdr, pi0, mu1, side){
  if(side == "one") {
    1 - pnorm(1/mu1*log((pi0/lfdr - pi0)/(1-pi0)) + mu1/2)
  } else if(side == "two") {
    tily = 2*pi0*(1-lfdr)/(lfdr*(1-pi0)*exp(-mu1^2/2))
    2*(1-pnorm(log((tily + sign(mu1)*sqrt(tily^2-4))/2)/mu1))
  }
}

lfdr.inverse_mvgauss_two <- function(lfdr){
  lfdr.inverse_mvgauss(lfdr, 0.9, 1, "two")
}
lfdr.inverse_mvgauss_two(0.94)


nonnull.cdf_mvgauss <- function(p, mu1, side) {
  if(side == "one") {
    pnorm(qnorm(1-p) - mu1,lower.tail = F)
  } else if(side == "two") {
    pnorm((qnorm(1-p/2) - mu1),lower.tail = F) + pnorm(qnorm(1-p/2))
  }
}

nonnull.cdf_mvgauss1 <- function(p){
  
}
  
  nonnull.cdf_mvgauss(p, 3, "one")



