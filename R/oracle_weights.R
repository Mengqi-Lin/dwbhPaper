
oracle.weights_mvgauss <- function(alpha, n_g, pi0_g, mu1_g, pi0Est = T, side) {
  if(alpha < 1e-10) {
    return(rep(0, sum(n_g)))
  }
  if(side == "right" | side == "left") {
    side <- "one"
  }
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
      R[i] = V[i] + n*(1-pi0)*nonnull.cdf_mvgauss(t, mu1, side)
    }
    return(sum(V)/sum(R))
  }
  #lambda_alpha <- uniroot(f= function(c){mFDR(c) - alpha}, lower = 1e-7, upper = 0.999)$root
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
  return(w)
}

lfdr.inverse_mvgauss <- function(lfdr, pi0, mu1, side){
  if(side == "one") {
    1 - pnorm(1/mu1*log((pi0/lfdr - pi0)/(1-pi0)) + mu1/2)
  } else if(side == "two") {
    tily = 2*pi0*(1-lfdr)/(lfdr*(1-pi0)*exp(-mu1^2/2))
    return(2*(1-pnorm(log((tily + sign(mu1)*sqrt(tily^2-4))/2)/mu1)))
  }
}

nonnull.cdf_mvgauss <- function(p, mu1, side) {
  if(side == "one") {
    pnorm(qnorm(1-p) - mu1,lower.tail = F)
  } else if(side == "two") {
    pnorm((qnorm(1-p/2) - mu1),lower.tail = F) + pnorm(-qnorm(1-p/2)-mu1)
  }
}



pval.ds_mvgauss <- function(p, mu1, pi0, side){
  if(side == "one") {
    pi0+(1-pi0)*(dnorm(qnorm(p, lower.tail = F)-mu1)/dnorm(qnorm(p, lower.tail = F)))
  } else if(side == "two") {
    pi0+(1-pi0)*((dnorm(qnorm(p/2, lower.tail = F)-mu1) + 
                    dnorm(qnorm(p/2, lower.tail = F)+mu1))
                 /2/dnorm(qnorm(p/2, lower.tail = F)))
  }
}


