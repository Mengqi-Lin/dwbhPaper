plot_results <- function(vals, methods, title,
                         cols, ltys, pchs,
                         ylim, ylab, lwd,
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
       cex = lwd,
       cex.axis = lwd,
       cex.lab = lwd,
       cex.main = lwd*1.2, 
       main = title, axes = FALSE)
  axis(side = 1, at = c(0, 0.05, 0.1, 0.15, 0.2))
  axis(side = 2)
  alpha_pt = 1:nalphas
  for (i in 1:length(methods)){
    points(alphalist, vals[i, ], lwd = lwd, cex = lwd, 
           type = 'l', col = cols[i], lty = ltys[i])
    points(alphalist, vals[i, alpha_pt], lwd = lwd, cex = lwd,
           col = cols[i], pch = pchs[i])
  }
  if (legend){
    legend("topleft", methods,
           col = cols, lty = ltys, pch = pchs, lwd = lwd,
           seg.len = 3, cex = cex.legend, bty = "n")
  }
}
