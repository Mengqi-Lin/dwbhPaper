#libraray()
source("dBH_utils.R")
source("utils.R")
genSigma <- function(n, rho = 0,
                     type = c("AR", "MA", "equi", "block", "iid"),
                     ifsqrt = F,
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
  if(!ifsqrt) {
    return(Sigma)
  }
  if (type == "iid") {
    sqrtSigma <- Sigma
  } else {
    sqrtSigma <- with(eigen(Sigma), vectors %*% (sqrt(values) * t(vectors)))
  }
  return(list(Sigma = Sigma, sqrtSigma = sqrtSigma))
}

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
  mu[inds] <- switch(mu_type,
                     `1` = rep(mu1, m),
                     `2` = rep(mu1, m) + rnorm(m, sd = 0.3), 
                     `3` = mu1 * (rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1)))
  return(mu)
}


## For each dwBH/dwBY, report wBH/wBY in the front. And also dwBH/dwBY init.
gen_methods_dwBH <- function(gamma,
                             weight_type,
                             skip_dBH2){
  expr_params <- expand.grid(
    gamma = gamma,
    weight_type = weight_type
  )
  
  # BH_methods <- sapply(weight_type, function(x){
  #   weight_type <- paste0("(", x, ")")
  #   tmp <- paste0("BH_", weight_type)
  #   c(tmp, paste0(tmp, "_safe"))
  # })
  # tmp <- "BH_(trivial)"
  # BH_methods <- c(tmp, paste0(tmp, "_safe"))
  # methods <- c(as.character(BH_methods))
  
  dBH_methods <- apply(expr_params, 1, function(x){
    if (is.na(x[1])){
      gamma <- "safe"
    } else {
      gamma <- x[1]
    }
    
    method1 <- ifelse(x[2]=="trivial", paste0("BH_", gamma), paste0("wBH_(", x[2],")_", gamma))
    method2 <- ifelse(x[2]=="trivial", paste0("dBH_", gamma), paste0("dwBH_(", x[2],")_", gamma))
    method3 <- ifelse(x[2]=="trivial", paste0("dBH_init_", gamma), paste0("dwBH_init_(", x[2],")_", gamma))
    c(method1, method2, method3)
  })
  methods <- as.character(dBH_methods)
  if (!skip_dBH2){
    dBH2_methods <- apply(expr_params, 1, function(x){
      if (is.na(x[1])){
        gamma <- "safe"
      } else {
        gamma <- x[1]
      }
      weight_type <- paste0("weighting(", x[2], ")")
      
      method1 <- paste0("dwBH2_", weight_type,
                        "_", gamma)
      method2 <- paste0("dwBH2_init_", weight_type,
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
                     sqrtSigma = NULL,
                     side, 
                     nreps){
  n <- sum(n_g)
  ngroups <- length(n_g)
  
  if (!(length(n_g) == length(mu1_g) & length(n_g) == length(pi1_g))){
    stop("Each group should have its corresponding pi1 and mu1.")
  }
  
  if(is.null(sqrtSigma)) {
    if(Sigma_type == "iid") {
      sqrtSigma <- diag(n)
    } else {
      Sigma <- genSigma(n, rho, Sigma_type)
      eigSigma <- eigen(Sigma)
      sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
    }
  } 
  
  mu <- c()
  groups <- c()
  for (j in 1: ngroups) {
    mu <- c(mu, genmu(n_g[j], pi1_g[j], mu1_g[j]))
    groups <- c(groups, rep(j , n_g[j]))
  }
  H0 <- mu == 0
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
  return(list(zvals = zvals, groups = groups, H0 = H0, sqrtSigma = sqrtSigma))
}


dwBH_mvgauss_groups_expr <- function(
  n_g, 
  mu1_g, 
  pi1_g, 
  mu_type = 2,
  posit_type = "fix",
  rho, 
  Sigma_type,
  side,
  alphas, 
  nreps, 
  weight_type,
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
    mu <- c(mu, genmu(n_g[j], pi1 = pi1_g[j], mu1 = mu1_g[j], mu_type = mu_type, posit_type = posit_type))
    groups <- c(groups, rep(j , n_g[j]))
  }
  if (side == "right"){
    mu <- abs(mu)
  } else if (side == "left"){
    mu <- -abs(mu)
  }
  
  H0 <- mu == 0
  nalphas <- length(alphas)
  Sigma <- genSigma(n, rho, Sigma_type)
  eigSigma <- eigen(Sigma)
  sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))
  
  methods <- gen_methods_dwBH(gamma, weight_type, skip_dBH2 = skip_dBH2)
  
  expr_params <- expand.grid(
    gamma = gamma,
    weight_type = weight_type
  )
  
  results <- lapply(1:nalphas, function(k){
    tmp <- matrix(NA, length(methods), nreps)
    rownames(tmp) <- methods
    return(list(alpha = alphas[k],
                FDP = tmp,
                power = tmp,
                secBH = tmp))
  })
  
  
  optimal_weights <- lapply(1:nalphas, function(k){
    lapply(1:nrow(expr_params), function(j){
      gam <- expr_params[j, "gamma"]
      gam <- ifelse(is.na(gam), 1/normalize(1:n), gam)
      oracle.weights_mvgauss(alpha = alphas[k]*gam, n_g = n_g, pi0_g = 1-pi1_g, mu1_g = mu1_g, pi0Est = T, side = side)
    })
  })
  
  pi1 <- sum(n_g*pi1_g)/n
  GBH_weights <- rep(pi1_g/(1-pi1_g)/pi1, n_g)
  trivial_weights <- rep(1, n)
  pb <- txtProgressBar(style=3)
  for (i in 1:nreps){
    zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
    pvals <- zvals_pvals(zvals, side)
    for (k in 1:nalphas){
      obj <- list()
      alpha <- alphas[k]
      
      ## BH rejections
      # avals <- 1:n
      # rejs_BH <- BH(pvals, alpha, avals, FALSE)
      # rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
      # obj <- c(obj, list(rejs_BH, rejs_BH_safe))
      
      # ## BC rejections
      # rejs_BC <- BC(pvals, alpha)
      # obj <- c(obj, list(rejs_BC))
      
      # ## Number of methods so far
      # nBHBC <- length(obj)
      
      ## dwBH rejections and wBH rejections
      
      for (j in 1:nrow(expr_params)){
        fac <- expr_params[j, 1]
        weight_type <- expr_params[j, 2]
        
        avals_type <- "BH"
        #avals <- 1:n
        #qvals <- qvals_BH_reshape(pvals, avals)
        if (is.na(fac)){
          gamma <- NULL
        } else {
          gamma <- fac
        }
        
        if(weight_type == "trivial") {
          weights <- trivial_weights
        } else if (weight_type == "GBH") {
          weights <- GBH_weights
        } else if (weight_type == "optimal") {
          weights <- optimal_weights[[k]][[j]]
        }
        
        rejs_dBH <- dBH_mvgauss(
          zvals = zvals,
          Sigma = Sigma,
          side = side,
          alpha = alpha,
          gamma = gamma, 
          niter = 1,
          tautype = "QC",
          weights = weights, 
          avals_type = avals_type)
        # rejs_dBH$maxq <- ifelse(
        #   length(rejs_dBH$initrejs) == 0, NA,
        #   max(qvals[rejs_dBH$initrejs] / alpha))
        rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
        rejs_BH <- list(rejs = rejs_dBH$Rhatrejs)
        obj <- c(obj, list(rejs_BH, rejs_dBH, rejs_dBH_init))
      }
      lengthdBH1 <- length(obj)
      if (!skip_dBH2){
        ## dBH2 rejections
        for (j in 1:nrow(expr_params)){
          fac <- expr_params[j, 1]
          weight_type  <- expr_params[j, 2]
          avals_type <- "BH"
          avals <- 1:n
          
          #qvals <- qvals_BH_reshape(pvals, avals)
          if (is.na(fac)){
            gamma <- NULL
          } else {
            gamma <- fac
          }
          
          if(weight_type == "trivial") {
            weights <- trivial_weights
          } else if (weight_type == "GBH") {
            weights <-GBH_weights
          } else if (weight_type == "optimal") {
            weights <- optimal_weights[[k]][[j]]
          }
          
          rejs_dBH2 <- dBH_mvgauss(
            zvals = zvals,
            Sigma = Sigma,
            side = side,
            alpha = alpha,
            gamma = gamma, 
            niter = 2,
            tautype = "QC",
            weights = weights,
            avals_type = avals_type,
            ...)
          # rejs_dBH2$maxq <- ifelse(
          #   length(rejs_dBH2$initrejs) == 0, NA,
          #   max(qvals[rejs_dBH2$initrejs] / alpha))
          rejs_dBH2_init <- list(rejs = rejs_dBH2$initrejs)
          obj <- c(obj, list(rejs_dBH2, rejs_dBH2_init))
        }
      }
      
      res <- sapply(obj, function(output){
        FDPpower(output$rejs, H0)
      })
      results[[k]]$FDP[, i] <- as.numeric(res[1, ])
      results[[k]]$power[, i] <- as.numeric(res[2, ])
      inds1 <- 1 + seq(1, lengthdBH1, 3)
      if(skip_dBH2) {
        inds <- inds1
      } else {
        inds2 <- seq(lengthdBH1+1, length(obj), 2)
        inds <- c(inds1, inds2)
      }
      results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
        output$secBH
      })
      # results[[k]]$sumWeights[inds, i] <- sapply(obj[inds], function(output){
      #   output$sumWeights
      # })
      # results[[k]]$qcap[inds, i] <- sapply(obj[inds], function(output){
      #   output$maxq
      # })
      setTxtProgressBar(pb, ((i-1)*nalphas + k)/(nreps*nalphas))
    }
  }
  
  close(pb)
  return(results)
}



postprocess <- function(res){
  summaryres <- lapply(res, function(re){
    FDR <- as.numeric(rowMeans(re$FDP))
    FDR <- round(FDR, 4)
    power <- as.numeric(rowMeans(re$power))
    secBH <- as.numeric(rowMeans(re$secBH))
    # qmax <- as.numeric(apply(re$qcap, 1, function(x){
    #   max(x, na.rm = TRUE)
    # }))
    # q99 <- as.numeric(apply(re$qcap, 1, function(x){
    #   quantile(x, 0.99, na.rm = TRUE)
    # }))
    # q95 <- as.numeric(apply(re$qcap, 1, function(x){
    #   quantile(x, 0.95, na.rm = TRUE)
    # }))
    methods <- rownames(re$power)
    df <- data.frame(method = methods,
                     FDR = FDR,
                     power = power,
                     secBH = secBH)
    # df_q <- data.frame(qmax = qmax,
    #                    q99 = q99,
    #                    q95 = q95)
    inds <- grep("^BH", methods)
    # inds2 <- grep("^BC", methods)
    # inds3 <- grep("^Knockoff", methods)
    #inds <- c(inds1, inds2, inds3)
    df1 <- df[inds, ]
    # df1_q <- df_q[inds, ]
    df1[, 5:6] <- NA
    names(df1)[5:6] <- c("FDR (init)", "power (init)")
    df2 <- df[-inds, ]
    m <- nrow(df2)
    df2_1 <- df2[seq(1, m, 2), ]
    df2_2 <- df2[seq(2, m, 2), ][, 2:3]
    names(df2_2) <- c("FDR (init)", "power (init)")
    df2 <- cbind(df2_1, df2_2)
    # df2_q <- df_q[-inds, ]
    # df2_q <- df2_q[seq(1, m, 2), ]
    
    df <- rbind(df1, df2)
    df <- df[, c(1, 2, 5, 3, 6, 4)]
    df$alpha <- re$alpha
    # df_q <- rbind(df1_q, df2_q)
    return(df)
  })
  do.call(rbind, summaryres)
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
