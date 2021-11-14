## fix mu_1
set.seed(1)
wBH_oracle_fix_mu1_twoside <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = rep(2, ng), 
                          pi1_g = seq(0.01, 0.1, length.out = ng), 
                          Sigma_type = "iid",
                          side = "two",
                          alphas = alphas, 
                          nreps = 1000, 
                          weight_type = weight_type,
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_fix_mu1_postres_twoside <- lapply(wBH_oracle_fix_mu1_twoside, function(x){
  wBH_postprocess(x)
})


FDR.filename <- "../figs/wBH_oracle_fixmu1_twoside_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_fix_mu1_postres_twoside[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_fixmu1_twoside_power.pdf"

pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_fix_mu1_postres_twoside[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.5), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()


## fix pi1
set.seed(1)
wBH_oracle_fixpi1_twoside <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 4, length.out = ng), 
                          pi1_g = rep(0.01, ng), 
                          Sigma_type = "iid",
                          side = "two",
                          alphas = alphas, 
                          nreps = 1000, 
                          weight_type = weight_type,
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_fixpi1_postres_twoside <- lapply(wBH_oracle_fixpi1_twoside, function(x){
  wBH_postprocess(x)
})


FDR.filename <- "../figs/wBH_oracle_fixpi1_twoside_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
titles <- c(paste0(ngroups[1], " groups "), paste0(ngroups[2], " groups "),
            paste0(ngroups[3], " groups "))
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_fixpi1_postres_twoside[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_fixpi1_twoside_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_fixpi1_postres_twoside[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.6), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()


# pi1 inc in mu1
set.seed(1)
wBH_oracle_inr_twoside <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 4, length.out = ng), 
                          pi1_g = seq(0.01, 0.05, length.out = ng), 
                          Sigma_type = "iid",
                          side = "two",
                          alphas = alphas, 
                          nreps = 1000, 
                          weight_type = weight_type,
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_inr_postres_twoside <- lapply(wBH_oracle_inr_twoside, function(x){
  wBH_postprocess(x)
})

FDR.filename <- "../figs/wBH_oracle_inc_twoside_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_inr_postres_twoside[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_inc_twoside_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_inr_postres_twoside[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()



# pi1 dec in mu1
set.seed(1)
wBH_oracle_dec_twoside <- lapply(ngroups, function(ng){
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = rev(seq(1, 4, length.out = ng)), 
                          pi1_g = seq(0.05, 0.1, length.out = ng), 
                          Sigma_type = "iid",
                          side = "two",
                          alphas = alphas, 
                          nreps = 1000, 
                          weight_type = weight_type,
                          skip_BY = F,
                          pi0Est = T)
})
wBH_oracle_dec_postres_twoside <- lapply(wBH_oracle_dec_twoside, function(x){
  wBH_postprocess(x)
})

FDR.filename <- "../figs/wBH_oracle_dec_twoside_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_dec_postres_twoside[[k]]$FDR
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = (k == 1), cex.legend = 1.1,
               alphas = alphas)
}
dev.off()

power.filename <- "../figs/wBH_oracle_dec_twoside_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_dec_postres_twoside[[k]]$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.5), ylab = "power",
               alphas = alphas,
               legend = F)
}
dev.off()



## tempting case
set.seed(1)
wBH_oracle_spe <-   wBH_mvgauss_groups_expr(n_g = c(1000, 1000, 1000, 1000), 
                                                   mu1_g = c(1, 3, 4, 2.5), 
                                                   pi1_g = c(0.3, 0.1, 0.01, 0.001), 
                                                   Sigma_type = "iid",
                                                   side = "two",
                                                   alphas = alphas, 
                                                   nreps = 100, 
                                                   weight_type = weight_type,
                                                   skip_BY = F,
                                                   pi0Est = T)

wBH_oracle_postres_spe <- wBH_postprocess(wBH_oracle_spe)

par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:1){
  power <- wBH_oracle_postres_spe$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1), ylab = "power",
               alphas = alphas,
               legend = F)
}





## tempting case??????
set.seed(1)
wBH_oracle_fixpi1_spe2 <-   wBH_mvgauss_groups_expr(n_g = c(1000, 1000), 
                                                   mu1_g = c(4, 1), 
                                                   pi1_g = c(0.01, 0.05), 
                                                   Sigma_type = "iid",
                                                   side = "two",
                                                   alphas = alphas, 
                                                   nreps = 100, 
                                                   weight_type = weight_type,
                                                   skip_BY = F,
                                                   pi0Est = T)

wBH_oracle_fixpi1_postres_spe2 <- wBH_postprocess(wBH_oracle_fixpi1_spe2)

par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:1){
  power <- wBH_oracle_fixpi1_postres_spe2$power
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 0.3), ylab = "power",
               alphas = alphas,
               legend = F)
}