source("wBH_utils.R")
wBH_mvgauss_groups_expr <- function(n_g, 
                                    mu1_g, 
                                    pi1_g, 
                                    rho, 
                                    Sigma_type,
                                    side, 
                                    skip_BY,
                                    alphas = c(0.05, 0.1, 0.15, 0.2), 
                                    nreps, 
                                    weight_type, 
                                    pi0Est){
  if (!(length(n_g) == length(mu1_g) & length(n_g) == length(pi1_g))){
    stop("Each group should have its corresponding pi1 and mu1.")
  }
  n <- sum(n_g)
  ngroups <- length(n_g)
  pi1 <- sum(n_g*pi1_g)/n
  pi0 <- 1-pi1
  
  # Trivial weights
  trivial_weights <- ifelse(pi0Est, rep(1/pi0, n), rep(1, n))
  
  ## GBH weights (can't correct without pi0)
  ## also generate mu and groups
  GBH_weights <- c()
  groups <- c()
  mu <- c()
  for (j in 1: ngroups) {
    mu <- c(mu, genmu(n_g[j], pi1 = pi1_g[j], mu1 = mu1_g[j], mu_type = 2))
    groups <- c(groups, rep(j , n_g[j]))
    GBH_weights <- c(GBH_weights, rep(pi1_g[j]/(1-pi1_g[j])/pi1, n_g[j]))
  }
  
  
  H0 <- mu == 0
  nalphas <- length(alphas)
  if(Sigma_type != "iid") {
    Sigma <- genSigma(n, rho, Sigma_type)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
  }
  
  ## optimal_weights for BH and BY
  optimal_weights <- list()
  optimal_weights_BY <- list()
  for(i in 1:nalphas) {
    alpha <- alphas[i]
    optimal_weights[[i]] <- oracle.weights_mvgauss(alpha = alpha, 
                                                   n_g = n_g, 
                                                   pi0_g = 1-pi1_g, 
                                                   mu1_g = mu1_g, 
                                                   pi0Est = pi0Est,
                                                   side = side)
    if(!skip_BY) {
      optimal_weights_BY[[i]] <- oracle.weights_mvgauss(alpha = alpha / normalize(1:n), 
                                                        n_g = n_g, 
                                                        pi0_g = 1-pi1_g, 
                                                        mu1_g = mu1_g, 
                                                        pi0Est =  pi0Est,
                                                        side = side)
    }
  }
  
  # gen methods
  expr_params <- expand.grid(
    weight_type = weight_type
  )
  
  wBH_methods <- apply(expr_params, 1, function(x){
    weight_type <- paste0("(", x[1], ")weighting")
    method1 <- paste0("wBH_", weight_type)
    if(!skip_BY) {
      method2 <- paste0("wBY_", weight_type)
    } else{
      method2 <- NULL
    }
    c(method1, method2)
  })
  
  methods <- c(as.character(wBH_methods))

  results <- lapply(1:nalphas, function(k){
    tmp <- matrix(NA, length(methods), nreps)
    rownames(tmp) <- methods
    return(list(alpha = alphas[k],
                FDP = tmp,
                power = tmp
                #,weights = tmp
                ))
  })
  
  pb <- txtProgressBar(style=3)
  for (i in 1:nreps){
    if(Sigma_type == "iid") {
      zvals <- mu + rnorm(n)
    } else{
      zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
    }
    pvals <- zvals_pvals(zvals, side)
    for (k in 1:nalphas){
      obj <- list()
      alpha <- alphas[k]
      for (j in 1:nrow(expr_params)){
        weight_type <- expr_params[j, 1]
        if(weight_type == "GBH" | weight_type == 2) {
          weights <- GBH_weights
        } else if(weight_type == "optimal" | weight_type == 3) {
          weights <- optimal_weights[[k]]
        } else if (weight_type == "trivial" | weight_type == 1) {
          weights <- trivial_weights
        } else if (weight_type == "adaptive_optimal_max" | weight_type == 4) {
          weights <- adaptive_optimal.weights(groups = groups, 
                                              zvals = zvals, 
                                              alpha = alpha, 
                                              side = side, 
                                              type = "max",
                                              pi0Est = pi0Est)
        } else if (weight_type == "adaptive_optimal_linexp" | weight_type == 5) {
          weights <- adaptive_optimal.weights(groups = groups, 
                                              zvals = zvals, 
                                              alpha = alpha, 
                                              side = side, 
                                              type = "lin_exp", 
                                              pi0Est = pi0Est)
        } else {
          stop("weight_type not appropriate")
        }
        rejs_BH <- BH(pvals = pvals/weights, alpha = alpha, reshape = FALSE)
        
        if(!skip_BY) {
          if(weight_type == "GBH" | weight_type == 2 | weight_type == "trivial" | weight_type == 1) {
            weights_BY <- weights
          } else if(weight_type == "optimal" | weight_type == 3) {
            weights_BY <- optimal_weights_BY[[k]]
          } else if (weight_type == "adaptive_optimal_max" | weight_type == 4) {
            weights_BY <- adaptive_optimal.weights(groups = groups, 
                                                   zvals = zvals, 
                                                   alpha = alpha/normalize(1:n), 
                                                   side = side,
                                                   pi0Est = pi0Est)
          } else if (weight_type == "adaptive_optimal_linexp" | weight_type == 5) {
            weights_BY <- adaptive_optimal.weights(groups = groups, 
                                                   zvals = zvals, 
                                                   alpha = alpha/normalize(1:n), 
                                                   side = side, 
                                                   type = "lin_exp",
                                                   pi0Est = pi0Est)
          } else {
            stop("weight_type not appropriate")
          }
          rejs_BH_safe <- BH(pvals = pvals/weights_BY, alpha = alpha, reshape = TRUE)
          obj <- c(obj, list(rejs_BH, rejs_BH_safe))
        } else {
          obj <- c(obj, list(rejs_BH))
        }
      }
      
      
      res <- sapply(obj, function(output){
        FDPpower(output$rejs, H0)
      })
      results[[k]]$FDP[, i] <- as.numeric(res[1, ])
      results[[k]]$power[, i] <- as.numeric(res[2, ])
      setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
    }
  }
  
  close(pb)
  return(results)
}

wBH_postprocess <- function(res){
  summaryres <- lapply(res, function(re){
    FDR <- as.numeric(rowMeans(re$FDP))
    FDR <- round(FDR, 4)
    power <- as.numeric(rowMeans(re$power))
    methods <- rownames(re$power)
    df <- data.frame(method = methods,
                     FDR = FDR,
                     power = power)
    df$alpha <- re$alpha
    return(df)
  })
  
  nmethods <- length(methods)
  nalphas <- length(methods)
  FDR <- matrix(NA, nmethods, nalphas)
  power <- matrix(NA, nmethods, nalphas)
  
  df <- do.call(cbind, summaryres)
  FDR <- df[, which(colnames(df)=="FDR")]
  power <- df[, which(colnames(df)=="power")]
  
  # methods <- rownames(df[,"method"])
  # rownames(FDR) <- methods
  # rownames(power) <- methods
  return(list(FDR = FDR, power = power))
}

