mFDR <- function(c, n_g, pi0_g, mu1_g) {
  V <- c()
  R <- c()
  for(i in 1:length(pi0_g)) {
    pi0 = pi0_g[i]
    mu1 = mu1_g[i]
    n = n_g[i]
    t = lfdr.inverse(c, pi0, mu1)
    V[i] = n*pi0*t
    R[i] = V[i] + n*(1-pi0)*nonnull.cdf(t, mu1)
  }
  return(sum(V)/sum(R))
}


pval.ds <- function(p, mu1, pi0){
  pi0+(1-pi0)*(dnorm(qnorm(p, lower.tail = F)-mu1)/dnorm(qnorm(p, lower.tail = F)))
}


pval.ds_mvgauss <- function(p, mu1, pi0, side){
  if(side == "one") {
    pi0+(1-pi0)*(dnorm(qnorm(p, lower.tail = F)-mu1)/dnorm(qnorm(p, lower.tail = F)))
  } else if(side == "two") {
    pi0+(1-pi0)*((dnorm(qnorm(p/2, lower.tail = F)-mu1) + dnorm(qnorm(p/2, lower.tail = F)+mu1))
                 /2/dnorm(qnorm(p/2, lower.tail = F)))
  }
}


lfdr_mvt <- function(p, pi0, mu1) {
  f1 <- dt(pt(p, df = 15, lower.tail = F) - mu1, df = 15)/dt(pt(p, df = 15, lower.tail = F), df = 15)
  return(pi0/(pi0+(1-pi0)*(f1)))
}
