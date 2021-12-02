delta = 5
sigma2 = 1

set.seed(0)
y = c(rep(0, 50), rep(2, 50), rep(0, 50), rep(-2, 50)) + rnorm(200, mean = 0, sd = sqrt(sigma2))
n = length(y)
gamma_set = c(0.01, 0.5, 1, 5, 10, 50)
DP_result = CV.search.DP.univar(y, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_DP_hat = unlist(DP_result$cpt_hat[[min_idx]])
cpt_DP_LR = local.refine.univar(cpt_DP_hat, y)

BS_object = BS.univar(y, 1, n, delta)
BS_result = thresholdBS(BS_object, tau = 2)
tree_BS = BS_result$BS_tree_trimmed
BS_result
cpt_BS_hat = sort(BS_result$cpt_hat[,1])
cpt_BS_LR = local.refine.univar(cpt_BS_hat, y)

BS_tune_result = tuneBSunivar(BS_object, y)
BS_tune_LR = local.refine.univar(BS_tune_result$cpt, y)

intervals = WBS.intervals(M = 300, lower = 1, upper = n)
WBS_object = WBS.univar(y, 1, n, intervals$Alpha, intervals$Beta, delta)
WBS_result = thresholdBS(WBS_object, tau = 3)
WBS_result
cpt_WBS_hat = sort(WBS_result$cpt_hat[,1])
cpt_WBS_LR = local.refine.univar(cpt_WBS_hat, y)

WBS_tune_result = tuneBSunivar(WBS_object, y)
WBS_tune_LR = local.refine.univar(WBS_tune_result$cpt, y)
