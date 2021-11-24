source("~/Documents/GitHub/dbh/R/compute_knots_mvgauss.R")
source("~/Documents/GitHub/dbh/R/dBH_mvgauss.R")
source("~/Documents/GitHub/dbh/R/dBH_mvgauss_qc.R" )
source("~/Documents/GitHub/dbh/R/adapOptimal_weights.R")
source("~/Documents/GitHub/dbh/R/dBH_mvgauss_qc_optimal.R")
source("~/Documents/GitHub/dbh/R/dBH_utils.R")
library(Rcpp)
sourceCpp("~/Documents/GitHub/dbh/src/RBH_homotopy.cpp")
source("expr_functions.R")
source("plot_function.R")


n_g = c(5000, 500)
MC = 10
alphas = seq(0.01, 0.1, length.out = 10)
skip_dBH2 = T
gamma = 1
weight_type = c("trivial", "optimal")
nreps = 100
## Different combinations of mu1 and pi1
MU1 <- matrix(c(1, 2.5, 2, 2, 1, 3), 2)
Pi1 <- matrix(c(0.1, 0.1, 0.01, 0.2, 0.01, 0.2), 2)


set.seed(1)
res <- lapply(1:3, function(i){
  dwBH_mvgauss_groups_expr(
    n_g = n_g, 
    mu1_g = MU1[, i], 
    pi1_g = Pi1[, i], 
    rho = 0.8, 
    Sigma_type = "AR",
    side = "right",
    alphas = alphas, 
    nreps = nreps, 
    MC = MC,
    weight_type = weight_type,
    gamma = 1,
    tautype = "QC",
    skip_dBH2 = skip_dBH2
  )
})

postres <- lapply(res, postprocess)
## Save data
filename <- "../data/mvgauss_adapOptimal_AR0.8_oneside_MC10_gamma1.RData"
save(postres, file = filename)

postres <- lapply(res, wBH_postprocess)
methods <- gen_methods(gamma = gamma, MC = MC, weight_type = weight_type, skip_dBH2 = T)
init_index <- grep("init", methods)
methods <- methods[-init_index]
title <- bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i]))

## Plot results
FDR.filename <- "../figs/FDR_mvgauss_adapOptimal_AR0.8_oneside_MC10_gamma1.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
cols <- c('black', 'grey','green', 'blue')
ltys <- c(2,2,1,1)
pchs <- c(1,1,0,2) 

pdf(FDR.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
lapply(1:3, function(i){
  plot_results(postres[[i]]$FDR[-init_index,], methods, 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 1.2,
               ylim = c(0, 0.12), ylab = "FDR",
               legend = T, cex.legend = 3,
               alphas = alphas)
  abline(a = 0, b=1, col = "red")
})
abline(a = 0, b=1, col = "red")
dev.off()

power.filename <- "../figs/power_mvgauss_adapOptimal_AR0.8_oneside_MC10_gamma1.pdf"
pdf(power.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
lapply(1:3, function(i){
  plot_results(postres[[i]]$power, methods, 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 1.2,
               ylim = c(0, 1), ylab = "power",
               legend = F, cex.legend = 3,
               alphas = alphas)
})
dev.off()