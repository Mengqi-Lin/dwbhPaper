plot_results <- function(vals, methods, title,
                         cols, ltys, pchs,
                         ylim, ylab,
                         legend = TRUE,
                         cex.legend = 1.2,
                         alphas){
  ### Plot results for each method
  alphalist <- alphas
  nalphas <- length(alphas)
  plot(0:1, 0:1, type = 'n',
       xlim = c(0, max(alphas)), ylim = ylim,
       xlab = expression(paste('Target FDR level ',alpha)),
       ylab = ylab,
       main = title, axes = FALSE)
  axis(side = 1, at = c(0, 0.1, 0.2, 0.3))
  axis(side = 2)
  alpha_pt = 1:nalphas
  for (i in 1:length(methods)){
    points(alphalist, vals[i, ],
           type = 'l', col = cols[i], lty = ltys[i])
    points(alphalist, vals[i, alpha_pt],
           col = cols[i], pch = pchs[i])
  }
  if (legend){
    legend("topleft", methods,
           col = cols, lty = ltys, pch = pchs,
           seg.len = 3, cex = cex.legend, bty = "n")
  }
}

## Plots for simulation 1
methods <- gen_methods.wBH(c("trivial", "GBH", "optimal"), skip_BY = F)
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('black', 'black','orange', 'orange', 'blue', 'blue')
ltys <- c(2,2,2,2,1,1)
pchs <- c(3,3,1,1,2,2)   


#load("../data/simul1.RData")
FDR.filename <- "../figs/wBH_oracle_inc_FDR.pdf"
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
alphas = c(0.05, 0.1, 0.2)
titles <- c(paste0(ngroups[1], " groups "), paste0(ngroups[2], " groups "),
            paste0(ngroups[3], " groups "))
ylim <- c(0, 0.25)
for (k in 1:3){
  FDP <- wBH_oracle_inr_postres[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = T, cex.legend = 1.1,
               alphas = c(0.05, 0.1, 0.2))
}
dev.off()



power.filename <- "../figs/wBH_oracle_inc_power2.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_oracle_inr_postres[[k]]$power
  legend <- FALSE
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1.05), ylab = "power",
               alphas = c(0.05, 0.1, 0.2),
               legend = T)
}

dev.off()



#load("../data/simul1.RData")
titles <- c("2 groups ", "3 groups ",
            "5 groups ")
alphas = c(0.05, 0.1, 0.2)
FDR.filename <- "../figs/wBH_adaptive_inc_FDR.pdf"
ylim <- c(0, 0.25)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  FDP <- wBH_adaptive_inr_postres[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods, titles[k], c(cols), c(ltys), c(pchs),
               ylim = ylim, ylab = "FDR",
               legend = legend)
}
dev.off()


power.filename <- "../figs/wBH_adaptive_inc_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_adaptive_inr_postres[[k]]$power
  legend <- FALSE
  plot_results(power, methods, titles[k], c(cols, "purple", "grey"), c(5, 6, ltys), c(5,6,pchs),
               ylim = c(0, 1.05), ylab = "power",
               legend = legend)
}
dev.off()







titles <- c("2 groups ", "3 groups ",
            "5 groups ")
alphas = c(0.05, 0.1, 0.2)
FDR.filename <- "../figs/wBH_adaptive_dec_FDR.pdf"
ylim <- c(0, 0.25)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  FDP <- wBH_adaptive_dec_postres[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods, titles[k], c(cols), c(ltys), c(pchs),
               ylim = ylim, ylab = "FDR",
               legend = legend)
}
dev.off()


power.filename <- "../figs/wBH_adaptive_dec_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- wBH_adaptive_dec_postres[[k]]$power
  legend <- FALSE
  plot_results(power, methods, titles[k], c(cols, "purple", "grey"), c(5, 6, ltys), c(5,6,pchs),
               ylim = c(0, 1.05), ylab = "power",
               legend = legend)
}
dev.off()