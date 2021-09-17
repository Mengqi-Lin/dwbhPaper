#library("dbh")

source("dBH_utils.R")
source("utils.R")

genSigma <- function(n, rho = 0,
                     type = c("AR", "MA", "equi", "block", "iid"),
                     bsize = 10){
  type <- type[1]
  if (type == "AR"){
    Sigma <- rho^(abs(outer(1:n, 1:n, "-")))
  } else if (type == "MA"){
    Sigma <- diag(n)
    Sigma[cbind(1:(n-1), 2:n)] <- rho
    Sigma[cbind(2:n, 1:(n-1))] <- rho
  } else if (type == "equi"){
    Sigma <- matrix(rho, n, n)
    diag(Sigma) <- 1
  } else if (type == "block"){
    m <- floor(n / bsize)
    if (n > bsize * m){
      warning("n is not divisible by bsize")
    }
    Sigma <- diag(n)
    blockSigma <- matrix(rho, bsize, bsize)
    diag(blockSigma) <- 1
    for (i in 1:m){
      inds <- ((i - 1) * bsize + 1):(i * bsize)
      Sigma[inds, inds] <- blockSigma
    }
  } else if (type == "iid"){
    Sigma <- diag(n)
  }
  return(Sigma)
}

#genSigma(8, 0.8, "AR")

genmu <- function(n, pi1, mu1,
                  posit_type = c("random", "fix"),
                  mu_type = 1:3){
  m <- ceiling(n * pi1)
  posit_type <- posit_type[1]
  mu_type <- mu_type[1]
  if (posit_type == "random"){
    inds <- seq(1, n, floor(1 / pi1))[1:m]
  } else if (posit_type == "fix"){
    inds <- 1:m
  }
  mu <- rep(0, n)
  altmu <- switch(mu_type,
                  `1` = rep(1, m),
                  `2` = rnorm(m),
                  `3` = rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1))
  mu[inds] <- mu1 * altmu
  mu
}



gen_methods <- function(gamma,
                        weight_type,
                        MC_type,
                        skip_dBH2){
  expr_params <- expand.grid(
    gamma = gamma,
    weight_type = weight_type,
    MC_type = MC_type
  )
  
  # BH_methods <- sapply(weight_type, function(x){
  #   weight_type <- paste0("weighting(", x, ")")
  #   tmp <- paste0("BH_", weight_type)
  #   c(tmp, paste0(tmp, "_safe"))
  # })
  methods <- c()

  dBH_methods <- apply(expr_params, 1, function(x){
    if (is.na(x[1])){
      gamma <- "safe"
    } else {
      gamma <- x[1]
    }
    
    weight_type <- paste0("weighting(", x[2], ")")

    MC_type <- x[3]
    
    method1 <- paste0("dwBH_", weight_type,
                      "_", MC_type,
                      "_", gamma)
    method2 <- paste0("dwBH_init_", weight_type,
                      "_", MC_type,
                      "_", gamma)
    method3 <- paste0("wBH_", weight_type,
                      "_", MC_type,
                      "_", gamma)
    c(method1, method2, method3)
  })
  methods <- c(methods, as.character(dBH_methods))
  if (!skip_dBH2){
    dBH2_methods <- apply(expr_params, 1, function(x){
      if (is.na(x[1])){
        gamma <- "safe"
      } else {
        gamma <- x[1]
      }
      weight_type <- paste0("weighting(", x[2], ")")
      tautype <- x[3]
      
      method1 <- paste0("dwBH2_", weight_type,
                        "_", MC_type,
                        "_", gamma)
      method2 <- paste0("dwBH2_init_", weight_type,
                        "_", MC_type,
                        "_", gamma)
      c(method1, method2)
    })
    methods <- c(methods,
                 as.character(dBH2_methods))
  }
  return(methods)
}

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
    stop("not idd")
    # eigSigma <- eigen(Sigma)
    # sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
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


