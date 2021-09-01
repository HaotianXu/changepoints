delta = 10
sigma2 = 1

set.seed(0)
y = c(rep(0, 50), rep(2, 50), rep(0, 50), rep(-2, 50)) + rnorm(200, mean = 0, sd = sqrt(sigma2))
n = length(y)
gamma.set = c(0.01, 0.5, 1, 5, 10, 50)
DP_result = CV.search.DP.univar(y, gamma.set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_DP_hat = unlist(DP_result$cpt_hat[[min_idx]])
cpt_DPlr_hat = local.refine.univar(cpt_hat, y, w = 1/3)

BS_result = threshold.BS(BS.univar(y, 1, n, delta), 2)
tree_BS = BS_result$BS_tree_trimmed
BS_result$change_points
cpt_BS_hat = sort(BS_result$change_points[,1])
cpt_BSlr_hat = local.refine.univar(cpt_BS_hat, y, w = 1/3)


intervals = WBS.intervals(M = 300, lower = 1, upper = n)
WBS_result = threshold.BS(WBS.univar(y, 1, n, intervals$Alpha, intervals$Beta, delta), 3)
tree_WBS = WBS_result$BS_tree_trimmed
WBS_result$change_points
cpt_WBS_hat = sort(WBS_result$change_points[,1])
cpt_BSlr_hat = local.refine.univar(cpt_WBS_hat, y, w = 1/3)

