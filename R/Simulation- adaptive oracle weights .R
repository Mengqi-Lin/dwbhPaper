
## Adaptive experiments

set.seed(1)
ngroups <- c(2, 3, 5)
wBH_adaptive_inr <- lapply(ngroups, function(ng) {
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 3, length.out = ng), 
                          pi1_g = seq(0.01, 0.1, length.out = ng), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 100, 
                          skip_BY = T,
                          weight_type = c("trivial",
                                          "optimal", 
                                          "adaptive_optimal_max"),
                          pi0Est = F)  
})

wBH_adaptive_inr_postres <- lapply(wBH_adaptive_inr, function(x) {
  wBH_postprocess(x)
})

methods <- gen_methods.wBH(c("trivial",
                             "optimal", 
                             "adaptive_optimal_max"), skip_BY = T)



set.seed(1)
ngroups <- c(2, 3, 5)
wBH_adaptive_dec <- lapply(ngroups, function(ng) {
  wBH_mvgauss_groups_expr(n_g = rep(1000, ng), 
                          mu1_g = seq(1, 4, length.out = ng), 
                          pi1_g = rev(seq(0.01, 0.1, length.out = ng)), 
                          Sigma_type = "iid",
                          side = "one",
                          alphas = c(0.05, 0.1, 0.2), 
                          nreps = 100, 
                          skip_BY = T,
                          weight_type = c("trivial",
                                          "optimal", 
                                          "adaptive_optimal_max"),
                          pi0Est = F)  
})

wBH_adaptive_dec_postres <- lapply(wBH_adaptive_dec, function(x) {
  wBH_postprocess(x)
})
