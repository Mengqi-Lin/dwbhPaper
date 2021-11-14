source("dBH_utils.R")

pvals_mvgauss <- function(zvals, Sigma, side){
    zvals <- zvals / sqrt(diag(Sigma))
    if (side == "right"){
        side <- "one"
    } else if (side == "left"){
        side <- "one"
        zvals <- -zvals
    }
    zvals_pvals(zvals, side)
}

pvals_mvt <- function(tvals, Sigma, df, side){
    tvals <- tvals / sqrt(diag(Sigma))
    if (side == "right"){
        side <- "one"
    } else if (side == "left"){
        side <- "one"
        tvals <- -tvals
    }
    tvals_pvals(tvals, df, side)
}

FDPpower <- function(rejs, H0){
    if (length(rejs) == 0){
        return(c(0, 0))
    }
    nrejs <- length(rejs)
    fp <- length(intersect(rejs, which(H0)))
    tp <- nrejs - fp
    FDP <- fp / nrejs
    power <- tp / sum(!H0)
    return(c(FDP, power))
}

gen_fast_AR <- function(n, rho){
    x <- rep(NA, n)
    x[1] <- rnorm(1)
    for (i in 2:n){
        x[i] <- x[i - 1] * rho + rnorm(1) * sqrt(1 - rho^2)
    }
    return(x)
}


