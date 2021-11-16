dwBH_mvgauss_groups_expr <- function(n_g, mu1_g, pi1_g, 
                                     rho, Sigma_type,
                                     side,
                                     alphas, nreps, weight_type,
                                     MC,
                                     gamma = 0.9,
                                     tautype = "QC",
                                     skip_dBH2 = TRUE,
                                     ...){
  if (!(length(n_g) == length(mu1_g) & length(n_g) == length(pi1_g))){
    stop("Each group should have its corresponding pi1 and mu1.")
  }
  n <- sum(n_g)
  ngroups <- length(n_g)
  groups <- c()
  mu <- c()
  for (j in 1: ngroups) {
    mu <- c(mu, genmu(n_g[j], pi1 = pi1_g[j], mu1 = mu1_g[j]))
    groups <- c(groups, rep(j , n_g[j]))
  }
  # if (side == "right"){
  #   mu <- abs(mu)
  # } else if (side == "left"){
  #   mu <- -abs(mu)
  # }
  H0 <- mu == 0
  nalphas <- length(alphas)
  Sigma <- genSigma(n, rho, Sigma_type)
  eigSigma <- eigen(Sigma)
  sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
  
  methods <- gen_methods(gamma, weight_type, MC,
                         skip_dBH2 = T)
  
  expr_params <- expand.grid(
    gamma = gamma,
    weight_type = weight_type,
    MC = MC
  )
  
  results <- lapply(1:nalphas, function(k){
    tmp <- matrix(NA, length(methods), nreps)
    rownames(tmp) <- methods
    return(list(alpha = alphas[k],
                FDP = tmp,
                power = tmp,
                secBH = tmp))
  })
  
  pb <- txtProgressBar(style=3)
  for (i in 1:nreps){
    zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
    pvals <- zvals_pvals(zvals, side)
    for (k in 1:nalphas){
      obj <- list()
      alpha <- alphas[k]
      avals <- 1:n
      ## BH rejections

      rejs_BH <- BH(pvals, alpha, avals, FALSE)
      rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
      obj <- c(obj, list(rejs_BH, rejs_BH_safe))
      
      ## Number of methods so far
      nBH <- length(obj)
      
      for (j in 1:nrow(expr_params)){
        fac <- expr_params[j, 1]
        weight_type <- expr_params[j, 2]
        MC <- expr_params[j, 3]

        
        avals_type <- "BH"
        avals <- 1:n
        qvals <- qvals_BH_reshape(pvals, avals)
        if (is.na(fac)){
          gamma <- NULL
        } else {
          gamma <- fac
        }
        rejs_dBH <- dBH_mvgauss(
          zvals = zvals,
          Sigma = Sigma,
          side = side,
          alpha = alpha,
          gamma = gamma, 
          covariates = groups,
          niter = 1,
          tautype = "QC",
          weight_type = weight_type,
          MC = MC,
          avals_type = avals_type)
        # rejs_dBH$maxq <- ifelse(
        #   length(rejs_dBH$initrejs) == 0, NA,
        #   max(qvals[rejs_dBH$initrejs] / alpha))
        rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
        obj <- c(obj, list(rejs_dBH, rejs_dBH_init))
      }
      
      if (!skip_dBH2){
        ## dBH2 rejections
        for (j in 1:nrow(expr_params)){
          fac <- expr_params[j, 1]
          weight_type  <- expr_params[j, 2]
          type <- expr_params[j, 3]
          avals_type <- "BH"
          avals <- 1:n
          
          qvals <- qvals_BH_reshape(pvals, avals)
          if (is.na(fac)){
            gamma <- NULL
          } else {
            gamma <- fac
          }
          rejs_dBH2 <- dBH_mvgauss(
            zvals = zvals,
            Sigma = Sigma,
            side = side,
            alpha = alpha,
            gamma = gamma, 
            covariates = groups,
            niter = 2,
            tautype = type,
            weight_type = weight_type,
            avals_type = avals_type,
            ...)
          rejs_dBH2$maxq <- ifelse(
            length(rejs_dBH2$initrejs) == 0, NA,
            max(qvals[rejs_dBH2$initrejs] / alpha))
          rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
          obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
        }
      }
      
      res <- sapply(obj, function(output){
        FDPpower(output$rejs, H0)
      })
      results[[k]]$FDP[, i] <- as.numeric(res[1, ])
      results[[k]]$power[, i] <- as.numeric(res[2, ])
      inds <- seq(nBH + 1, length(methods), 2)
      results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
        output$secBH
      })
      # results[[k]]$qcap[inds, i] <- sapply(obj[inds], function(output){
      #   output$maxq
      # })
      setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
    }
  }
  
  close(pb)
  return(results)
}



