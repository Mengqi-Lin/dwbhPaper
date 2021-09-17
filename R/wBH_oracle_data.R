nGroups <- c(10, 20, 30)
result1 <- list()
set.seed(1)
for(i in 1:length(nGroups)) {
  ngroups <- nGroups[i]
  result1[[i]] <- wBH_mvgauss_groups_expr(n_g = rep(1000, ngroups), 
                                         mu1_g = rev(seq(1, 6, length.out = ngroups)), 
                                         pi1_g = seq(0.01, 0.3, length.out = ngroups), 
                                         Sigma_type = "iid",
                                         side = "one",
                                         alphas = c(0.05, 0.1, 0.15, 0.2), 
                                         nreps = 1000, 
                                         weight_type = c("trivial", "GBH", "optimal")) %>% 
                wBH_postprocess
}
filename <- "../data/wBH_dec.RData"
save(file = filename, result1)

set.seed(1)
for(i in 1:length(nGroups)) {
  ngroups <- nGroups[i]
  result[[i]] <- wBH_mvgauss_groups_expr(n_g = rep(1000, ngroups), 
                                         mu1_g = seq(1, 6, length.out = ngroups), 
                                         pi1_g = seq(0.01, 0.3, length.out = ngroups), 
                                         Sigma_type = "iid",
                                         side = "one",
                                         alphas = c(0.05, 0.1, 0.15, 0.2), 
                                         nreps = 100, 
                                         weight_type = c("trivial", "GBH", "optimal")) %>% 
                wBH_postprocess
}
filename <- "../data/wBH_inc.RData"
save(file = filename, result)






result3 <- list()
set.seed(1)
for(i in 1:length(nGroups)) {
  ngroups <- nGroups[i]
  result3[[i]] <- wBH_mvgauss_groups_expr(n_g = rep(1000, ngroups), 
                                          mu1_g = rev(seq(1, 6, length.out = ngroups)), 
                                          pi1_g = seq(0.01, 0.3, length.out = ngroups), 
                                          Sigma_type = "iid",
                                          side = "one",
                                          alphas = c(0.05, 0.1, 0.15, 0.2), 
                                          nreps = 100, 
                                          weight_type = c("trivial", "GBH", "optimal", "adaptive_optimal")) %>% 
    wBH_postprocess
}
filename <- "../data/wBH_dec_adaptive.RData"
save(file = filename, result3)

