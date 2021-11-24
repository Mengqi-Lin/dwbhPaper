#library("dbh")

source("dBH_utils.R")
source("utils.R")
source("oracle_weights.R")
gen_data <- function(n_g, mu1_g, pi1_g, 
                     rho, Sigma_type,
                     side, 
                     nreps){
  n <- sum(n_g)
  ngroups <- length(n_g)
  if (!(length(n_g) == length(mu1_g) & length(n_g) == length(pi1_g))){
    stop("Each group should have its corresponding pi1 and mu1.")
  }
  Sigma <- genSigma(n, rho, Sigma_type)
  if(Sigma_type == "iid") {
    sqrtSigma <- Sigma
  } else {
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
  }
  
  mu <- c()
  groups <- c()
  for (j in 1: ngroups) {
    mu <- c(mu, genmu(n_g[j], pi1_g[j], mu1_g[j]))
    groups <- c(groups, rep(j , n_g[j]))
  }
  H0 <- mu == 0
  #oracle.weights <- oracle.weights(alpha, n_g, pi0_g = 1- pi1_g, mu1_g = mu1_g)
  zvals <- list()
  for (i in 1:nreps){
    if (side == "right"){
      mu <- abs(mu)
    } else if (side == "left"){
      mu <- -abs(mu)
    }
    zvals[[i]] <- as.numeric(mu + sqrtSigma %*% rnorm(n))
    # zvals[[i]] <- zvals[[i]] / sqrt(diag(Sigma))
    # if (!is.null(Sigma) && any(diag(Sigma) != 1)){
    #   Sigma <- cov2cor(Sigma)
    # }
  }
  return(list(zvals = zvals, groups = groups, H0 = H0, Sigma = Sigma))
}

