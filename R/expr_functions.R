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
                        MC,
                        skip_dBH2){
    expr_params <- expand.grid(
        gamma = gamma,
        weight_type = weight_type,
        MC = MC
    )
    
    methods <- c("BH", "BY")
    
    dBH_methods <- apply(expr_params, 1, function(x){
        if (is.na(x[1])){
            gamma <- "safe"
        } else {
            gamma <- x[1]
        }
        
        weight_type <- paste0("weighting(", x[2], ")")
        
        MC <- x[3]
        
        method1 <- paste0("dwBH_", weight_type,
                          "_", MC,
                          "_", gamma)
        method2 <- paste0("dwBH_init_", weight_type,
                          "_", MC,
                          "_", gamma)
        c(method1, method2)
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
                              "_", MC,
                              "_", gamma)
            method2 <- paste0("dwBH2_init_", weight_type,
                              "_", MC,
                              "_", gamma)
            c(method1, method2)
        })
        methods <- c(methods,
                     as.character(dBH2_methods))
    }
    return(methods)
}

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

dBH_mvt_expr <- function(n, df, mu1, pi1,
                         mu_posit_type, mu_size_type,
                         rho, Sigma_type,
                         side,
                         alphas, nreps,
                         gamma = 0.9,
                         geom_fac = 2,
                         tautype = "QC",
                         skip_dBH2 = TRUE,
                         ...){
    nalphas <- length(alphas)
    Sigma <- genSigma(n, rho, Sigma_type)
    eigSigma <- eigen(Sigma)
    sqrtSigma <- with(eigSigma, vectors %*% (sqrt(values) * t(vectors)))

    methods <- gen_methods(gamma, geom_fac, tautype,
                           TRUE, skip_dBH2)
    expr_params <- expand.grid(
        gamma = gamma,
        geom_fac = geom_fac,
        tautype = tautype
    )

    results <- lapply(1:nalphas, function(k){
        tmp <- matrix(NA, length(methods), nreps)
        rownames(tmp) <- methods
        return(list(alpha = alphas[k],
                    FDP = tmp,
                    power = tmp,
                    secBH = tmp,
                    qcap = tmp))
    })
    
    pb <- txtProgressBar(style=3)
    for (i in 1:nreps){
        mu <- genmu(n, pi1, mu1, mu_posit_type, mu_size_type)
        if (side == "right"){
            mu <- abs(mu)
        } else if (side == "left"){
            mu <- -abs(mu)
        }
        H0 <- mu == 0
        zvals <- as.numeric(mu + sqrtSigma %*% rnorm(n))
        sigmahat <- sqrt(rchisq(1, df = df) / df)
        tvals <- zvals / sigmahat
        pvals <- pvals_mvt(tvals, Sigma, df, side)

        for (k in 1:nalphas){
            obj <- list()
            alpha <- alphas[k]            

            ## BH rejections
            for (x in union(NA, geom_fac)){
                if (is.na(x)){
                    avals <- 1:n
                } else {
                    avals <- geom_avals(x, n)
                }
                rejs_BH <- BH(pvals, alpha, avals, FALSE)
                rejs_BH_safe <- BH(pvals, alpha, avals, TRUE)
                obj <- c(obj, list(rejs_BH, rejs_BH_safe))
            }

            ## BC rejections
            rejs_BC <- BC(pvals, alpha)
            obj <- c(obj, list(rejs_BC))

            ## Number of methods so far
            nBHBC <- length(obj)

            ## dBH rejections
            for (j in 1:nrow(expr_params)){
                fac <- expr_params[j, 1]
                x <- expr_params[j, 2]
                type <- expr_params[j, 3]
                if (is.na(x)){
                    avals_type <- "BH"
                    avals <- 1:n                    
                } else {
                    avals_type <- "geom"
                    avals <- geom_avals(x, n)
                }
                qvals <- qvals_BH_reshape(pvals, avals)
                if (is.na(fac)){
                    gamma <- NULL
                } else {
                    gamma <- fac
                }
                rejs_dBH <- dBH_mvt(
                    tvals = tvals,
                    df = df,
                    Sigma = Sigma,
                    side = side,
                    alpha = alpha,
                    gamma = gamma, 
                    niter = 1,
                    tautype = type,
                    avals_type = avals_type,
                    geom_fac = x, ...)
                rejs_dBH$maxq <- ifelse(
                    length(rejs_dBH$initrejs) == 0, NA,
                    max(qvals[rejs_dBH$initrejs] / alpha))
                rejs_dBH_init <- list(rejs = rejs_dBH$initrejs)
                obj <- c(obj, list(rejs_dBH, rejs_dBH_init))
            }

            if (!skip_dBH2){
                ## dBH2 rejections
                for (j in 1:nrow(expr_params)){
                    fac <- expr_params[j, 1]
                    x <- expr_params[j, 2]
                    type <- expr_params[j, 3]
                    if (is.na(x)){
                        avals_type <- "BH"
                        avals <- 1:n                        
                    } else {
                        avals_type <- "geom"
                        avals <- geom_avals(x, n)
                    }
                    if (is.na(fac)){
                        gamma <- NULL
                    } else {
                        gamma <- fac
                    }
                    qvals <- qvals_BH_reshape(pvals, avals)
                    rejs_dBH2 <- dBH_mvt(
                        tvals = tvals,
                        df = df,
                        Sigma = Sigma,
                        side = side,
                        alpha = alpha,
                        gamma = gamma,
                        niter = 2,
                        tautype = type,
                        avals_type = avals_type,
                        geom_fac = x, ...)
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
            inds <- seq(nBHBC + 1, length(methods), 2)
            results[[k]]$secBH[inds, i] <- sapply(obj[inds], function(output){
                output$secBH
            })
            results[[k]]$qcap[inds, i] <- sapply(obj[inds], function(output){
                output$maxq
            })
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
        
        methods <- rownames(re$power)
        df <- data.frame(method = methods,
                         FDR = FDR,
                         power = power,
                         secBH = secBH)
        
        inds1 <- grep("^BH", methods)
        inds2 <- grep("^BY", methods)
        inds <- c(inds1, inds2)
        df1 <- df[inds, ]
        df1[, 5:6] <- NA
        names(df1)[5:6] <- c("FDR (init)", "power (init)")
        df2 <- df[-inds, ]
        m <- nrow(df2)
        df2_1 <- df2[seq(1, m, 2), ]
        df2_2 <- df2[seq(2, m, 2), ][, 2:3]
        names(df2_2) <- c("FDR (init)", "power (init)")
        df2 <- cbind(df2_1, df2_2)
        
        df <- rbind(df1, df2)
        df <- df[, c(1, 2, 5, 3, 6, 4)]
        df$alpha <- re$alpha
        return(cbind(df))
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

aggregate_expr <- function(objlist, alphas){
    res <- list()
    for (k in 1:length(alphas)){
        FDP <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$FDP
        }))
        power <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$power
        }))
        secBH <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$secBH
        }))
        qcap <- do.call(cbind, lapply(objlist, function(x){
            x[[k]]$qcap
        }))
        res[[k]] <- list(alpha = alphas[k], FDP = FDP, power = power, secBH = secBH, qcap = qcap)
    }
    return(res)
}
