dwBH_mvgauss_groups_expr <- function(n_g, mu1_g, pi1_g, 
                             mu_size_type,
                             rho, Sigma_type,
                             side,
                             alphas, nreps, weight_type,
                             MC_type,
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
    mu <- c(mu, genmu(n_g[j], pi1 = pi1_g[j], mu1 = mu1_g[j], posit_type = "fix"))
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
  
  methods <- gen_methods(gamma, weight_type, MC_type,
                         skip_dBH2 = T)
  expr_params <- expand.grid(
    gamma = gamma,
    weight_type = weight_type,
    MC_type = MC_type
  )
  
  results <- lapply(1:nalphas, function(k){
    tmp <- matrix(NA, length(methods), nreps)
    rownames(tmp) <- methods
    return(list(alpha = alphas[k],
                FDP = tmp,
                power = tmp,
                secBH = tmp,
                sumWeights = tmp))
  })
  
  pb <- txtProgressBar(style=3)
  for (i in 1:nreps){
    zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
    pvals <- zvals_pvals(zvals, side)
    for (k in 1:nalphas){
      obj <- list()
      alpha <- alphas[k]
      
      # ## BH rejections
      # for (x in union(NA, geom_fac)){
      #   if (is.na(x)){
      #     avals <- 1:n
      #   } else {
      #     avals <- geom_avals(x, n)
      #   }
      #   rejs_BH <- BH(pvals, alpha, avals, FALSE)
      #   rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
      #   obj <- c(obj, list(rejs_BH, rejs_BH_safe))
      # }
      # 
      # ## BC rejections
      # rejs_BC <- BC(pvals, alpha)
      # obj <- c(obj, list(rejs_BC))
      
      # ## Number of methods so far
      # nBHBC <- length(obj)
      
      ## dwBH rejections and wBH rejections
      for (j in 1:nrow(expr_params)){
        fac <- expr_params[j, 1]
        weight_type <- expr_params[j, 2]
        MC_type <- expr_params[j, 3]
        
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
          MC_type = MC_type,
          pi0_oracle = 1-pi1_g,
          mu1_oracle = mu1_g,
          avals_type = avals_type)
        # rejs_dBH$maxq <- ifelse(
        #   length(rejs_dBH$initrejs) == 0, NA,
        #   max(qvals[rejs_dBH$initrejs] / alpha))
        rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
        rejs_wBH <- list(rejs = rejs_dBH$dBH_rej0)
        obj <- c(obj, list(rejs_dBH, rejs_dBH_init, rejs_wBH))
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
            pi0_oracle = 1-pi1_g,
            mu1_oracle = mu1_g,
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
      inds <- seq(1, length(obj), 3)
      results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
        output$secBH
      })
      results[[k]]$sumWeights[inds, i] <- sapply(obj[inds], function(output){
        output$sumWeights
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



