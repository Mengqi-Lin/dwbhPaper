plot_results <- function(vals, methods, title,
                         cols, ltys, pchs,
                         ylim, ylab,
                         legend = TRUE,
                         cex.legend = 1.2){
  ### Plot results for each method
  alphalist <- alphas
  
  plot(0:1, 0:1, type = 'n',
       xlim = range(alphalist), ylim = ylim,
       xlab = expression(paste('Target FDR level ',alpha)),
       ylab = ylab,
       main = title, axes = FALSE)
  axis(side = 1, at = c(0, 0.1, 0.2, 0.3))
  axis(side = 2)
  alpha_pt = 1:nalphas
  for (i in 1:length(methods)){
    points(alphalist, vals[i, ],
           type = 'l', col = cols[i], lty = ltys[i])
    points(alphalist[alpha_pt], vals[i, alpha_pt],
           col = cols[i], pch = pchs[i])
  }
  if (legend){
    legend("topleft", methods,
           col = cols, lty = ltys, pch = pchs,
           seg.len = 3, cex = cex.legend, bty = "n")
  }
}


result <- list()
result[[1]] <- wBHpostres7
result[[2]] <- wBHpostres9
result[[3]] <- wBHpostres10

## Plots for simulation 1
methods <- unique(methods1)
## Blue for non-adaptive BH-type methods: BH, Storey, BC; Orange for IHW; Black for SABHA; Red for AdaPT
cols <- c('black', 'red','orange', 'green', 'blue', 'blue')
ltys <- c(1,2,3,2,3,1)
pchs <- c(2,3,1,1,2,1)   



#load("../data/simul1.RData")
titles <- c("10 groups ", "20 groups ",
            "30 groups ")
FDR.filename <- "../figs/wBH_oracle_inc_FDR.pdf"
ylim <- c(0, 0.25)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  FDP <- result[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = legend, cex.legend = 1.1)
}
dev.off()


power.filename <- "../figs/wBH_oracle_inc_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- result[[k]]$power
  legend <- FALSE
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1.05), ylab = "power",
               legend = legend)
}
dev.off()




#load("../data/simul1.RData")
titles <- c("10 groups ", "20 groups ",
            "30 groups ")
FDR.filename <- "../figs/wBH_oracle_dec_FDR.pdf"
ylim <- c(0, 0.25)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  FDP <- result1[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods, titles[k], cols, ltys, pchs,
               ylim = ylim, ylab = "FDR",
               legend = legend, cex.legend = 1.1)
}
dev.off()


power.filename <- "../figs/wBH_oracle_dec_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- result1[[k]]$power
  legend <- FALSE
  plot_results(power, methods, titles[k], cols, ltys, pchs,
               ylim = c(0, 1.05), ylab = "power",
               legend = legend)
}
dev.off()






#load("../data/simul1.RData")
titles <- c("10 groups ", "20 groups ",
            "30 groups ")
FDR.filename <- "../figs/wBH_adaptive_dec_FDR.pdf"
ylim <- c(0, 0.25)
pdf(FDR.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  FDP <- result3[[k]]$FDR
  legend <- (k == 1)
  plot_results(FDP, methods2, titles[k], c(cols, "purple", "grey"), c(5, 6, ltys), c(5,6,pchs),
               ylim = ylim, ylab = "FDR",
               legend = legend, cex.legend = 1.1)
}
dev.off()


power.filename <- "../figs/wBH_adaptive_dec_power.pdf"
pdf(power.filename, width = 9, height = 2.8)
par(mfrow = c(1, 3), mar = c(4, 5, 2, 2), oma=c(0, 0, 0, 5), cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.7)
for (k in 1:3){
  power <- result3[[k]]$power
  legend <- FALSE
  plot_results(power, methods2, titles[k], c(cols, "purple", "grey"), c(5, 6, ltys), c(5,6,pchs),
               ylim = c(0, 1.05), ylab = "power",
               legend = legend)
}
dev.off()

methods2 <- c(methods, "wBH_(AdapOpt)weighting", "wBY_(AdapOpt)weighting")
