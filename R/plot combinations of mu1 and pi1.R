source("utils.R")
source("dBH_utils.R")
source("wBH_utils.R")
source("oracle_weights.R")
source("wBH_mvgauss_groups_expr.R")




set.seed(1)
wBH_oracle_comb <- lapply(1:ncol(comb), function(i){
  wBH_mvgauss_groups_expr(n_g = rep(1000, 2), 
                          mu1_g = comb[3:4, i], 
                          pi1_g = comb[1:2, i], 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 1000, 
                          weight_type = c("trivial", "GBH", "optimal"),
                          skip_BY = F,
                          pi0Est = T)
})


wBH_oracle_comb_postres <- lapply(wBH_oracle_comb, function(x){
  wBH_postprocess(x)
})



FDR.filename <- "../figs/comb_FDR.pdf"
pdf(FDR.filename, width = 9, height = 6)
par(mfrow = c(2, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- paste0(comb[1,], " vs ", comb[2,], ", ", comb[3,], " vs ",comb[4,])
ylim <- c(0, 0.25)
for (k in 1:6){
  FDP <- wBH_oracle_comb_postres[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/comb_power.pdf"
pdf(power.filename, width = 9, height = 3)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:6){
  power <- wBH_oracle_comb_postres[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()
