source("utils.R")
source("dBH_utils.R")
source("oracle_weights.R")
source("dwBH_expr_functions.R")
source("wBH_oracle_plots.R")
library(EbayesThresh)

alphas = c(0, 0.05, 0.1, 0.15, 0.2)
skip_dBH2 = T

## Different combinations of mu1 and pi1
MU1 <- matrix(c(1, 2.5, 2, 2, 1, 3), 2)
Pi1 <- matrix(c(0.1, 0.1, 0.01, 0.2, 0.01, 0.2), 2)

set.seed(1)
res <- lapply(1:3, function(i){
  dwBH_mvgauss_groups_expr(
    n_g = c(1000, 1000), 
    mu1_g = MU1[, i], 
    pi1_g = Pi1[, i], 
    mu_type = 1,
    posit_type = "fix",
    rho = 0.8, 
    Sigma_type = "AR",
    side = "right",
    alphas = alphas, 
    nreps = 1000, 
    weight_type = c("trivial", "optimal"),
    gamma = c(1, NA),
    tautype = "QC",
    skip_dBH2 = skip_dBH2
  )
})

postres <- lapply(res, function(i){
  wBH_postprocess(i)
})
dwBH_posetres <- lapply(res, postprocess)

filename <- "../data/dwBH_mvgauss_AR0.8.RData"
save(res, file = filename)

cols <- c('black', 'black','orange', 'orange', 'blue', 'blue')
ltys <- c(1,2,1,2)
pchs <- c(1,2,1,2) 

methods <- gen_methods_dwBH(gamma = c(1, NA), weight_type = c("trivial", "optimal"), skip_dBH2 = T)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
BY_ind <- grep("safe", methods)
BH_ind <- grep(1, methods)
init_ind <- grep("init", methods)
BH_final_ind <- setdiff(BH_ind, init_ind)
BY_final_ind <- setdiff(BY_ind, init_ind)
## FDR plots
FDR.filename <- "../figs/dwBH_FDR_AR0.8_gaussian.pdf"
pdf(FDR.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
lapply(1:3, function(i){
  plot_results(postres[[i]]$FDR[BH_final_ind,], methods[BH_final_ind], 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 2.2,
               ylim = c(0, 0.25), ylab = "FDR",
               legend = T, cex.legend = 3,
               alphas = alphas)
  abline(a = 0, b=1, col = "red")
})
dev.off()

FDR.filename <- "../figs/dwBH_power_AR0.8_gaussian.pdf"
pdf(FDR.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
## power plots
lapply(1:3, function(i){
  plot_results(postres[[i]]$power[BH_final_ind,], 
               methods[BH_final_ind], 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 2.2,
               ylim = c(0,1), ylab = "power",
               legend = T, cex.legend = 1.1,
               alphas = alphas)
})
dev.off()


## FDR plots
FDR.filename <- "../figs/dwBY_FDR_AR0.8_gaussian.pdf"
pdf(FDR.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
lapply(1:3, function(i){
  plot_results(postres[[i]]$FDR[BY_final_ind,], methods[BY_final_ind], 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 2.2,
               ylim = c(0, 0.25), ylab ="FDR",
               legend = T, cex.legend = 3,
               alphas = alphas)
  abline(a=0, b= 1, col = "red")
})
dev.off()

FDR.filename <- "../figs/dwBY_power_AR0.8_gaussian.pdf"
pdf(FDR.filename, width = 21, height = 7)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2))
## power plots
lapply(1:3, function(i){
  plot_results(postres[[i]]$power[BY_final_ind,], 
               methods[BY_final_ind], 
               title = bquote(mu[1]~ ": " ~ .(MU1[1, i])~" vs "~ .(MU1[2, i])~", "~ pi[1]~": "~ .(Pi1[1, i])~" vs "~ .(Pi1[2, i])),
               cols = cols, ltys = ltys, pchs = pchs, lwd = 2.2,
               ylim = c(0,1), ylab = "power",
               legend = T, cex.legend = 3,
               alphas = alphas)
})
dev.off()


