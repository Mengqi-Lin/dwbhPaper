source("utils.R")
source("dBH_utils.R")
source("wBH_utils.R")
source("oracle_weights.R")
source("wBH_mvgauss_groups_expr.R")

ngroups <- c(2, 5, 10)
alphas <-  c(0.05, 0.1, 0.2)
weight_type = c("trivial", "GBH", "optimal")
skip_BY = F
methods <- gen_methods.wBH(c("trivial", "GBH", "optimal"), skip_BY = skip_BY)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- paste0(ngroups, " groups ")

## fix mu_1
set.seed(1)
wBH_oracle_fix_mu1 <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = rep(2, ng), 
                          pi1_g = seq(0.01, 0.1, length.out = ng), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 1000, 
                          weight_type = c("trivial", "GBH", "optimal"),
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_fix_mu1_postres <- lapply(wBH_oracle_fix_mu1, function(x){
  wBH_postprocess(x)
})


FDR.filename <- "../figs/wBH_oracle_fixmu1_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- c(paste0(ngroups[1], " groups "), paste0(ngroups[2], " groups "),
            paste0(ngroups[3], " groups "))
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_fix_mu1_postres[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_fixmu1_power.pdf"

pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_fix_mu1_postres[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.5), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()


## fix pi1
set.seed(1)
wBH_oracle_fixpi1 <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 4, length.out = ng), 
                          pi1_g = rep(0.01, ng), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 1000, 
                          weight_type = c("trivial", "GBH", "optimal"),
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_fixpi1_postres <- lapply(wBH_oracle_fixpi1, function(x){
  wBH_postprocess(x)
})


FDR.filename <- "../figs/wBH_oracle_fixpi1_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- c(paste0(ngroups[1], " groups "), paste0(ngroups[2], " groups "),
            paste0(ngroups[3], " groups "))
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_fixpi1_postres[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_fixpi1_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_fixpi1_postres[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.6), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()


# pi1 inc in mu1
set.seed(1)
ngroups <- c(2, 5, 10)
wBH_oracle_inr <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 4, length.out = ng), 
                          pi1_g = seq(0.01, 0.05, length.out = ng), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 1000, 
                          weight_type = c("trivial", "GBH", "optimal"),
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_inr_postres <- lapply(wBH_oracle_inr, function(x){
  wBH_postprocess(x)
})

FDR.filename <- "../figs/wBH_oracle_inc_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- c(paste0(ngroups[1], " groups "), paste0(ngroups[2], " groups "),
            paste0(ngroups[3], " groups "))
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_inr_postres[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_inc_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_inr_postres[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()



# pi1 dec in mu1
set.seed(1)
ngroups <- c(2, 5, 10)
wBH_oracle_dec <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = rev(seq(1, 4, length.out = ng)), 
                          pi1_g = seq(0.01, 0.05, length.out = ng), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 1000, 
                          weight_type = c("trivial", "GBH", "optimal"),
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_dec_postres <- lapply(wBH_oracle_dec, function(x){
  wBH_postprocess(x)
})

FDR.filename <- "../figs/wBH_oracle_dec_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_dec_postres[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_dec_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_dec_postres[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.5), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()
