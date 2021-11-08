d0 = 5# the number of nonzero elements = p-sparsity
RR = 5## replication
sigma = 1## variance of data
kappa = 5## spectral norm
delta = 3 ## minimal spacing to help DP more robust
n = 200## sample size, need to be the multiple of 2*gridsize
p = 50 ## dimensionality
cpt_true = c(80, 170)
gamma_dp_set = c(0.01,0.5,1,5,10,50)
lambda_dp_set = c(0.01,0.1,1,1.5)#c(0.01,0.1,1)#c(0.01,0.1,1)
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  data = simu.change.regression(d0, cpt_true, p, n, sigma, kappa)
  X = data$X
  y = data$y
  cpt_true = data$cpt_true
  result = CV.search.DP.regression(y, X, gamma_dp_set, lambda_dp_set, delta, eps = 0.001)
  ##output a table which contains estimated change points, number of estimated change points, 
  ##validation loss and training error.
  min_idx = as.vector(arrayInd(which.min(result$test_error), dim(result$test_error)))
  cpt_est_init[[time]] = unlist(result$cpt_hat[min_idx[1], min_idx[2]])
  error_init[time] = Hausdorff.dist(cpt_est_init[[time]], cpt_true)
  len_est_init[time] = unlist(result$K_hat[min_idx[1], min_idx[2]])
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
cpt_est_init ## estimated change points
error_init ##Hausdorff
len_est_init ## number of estimated change points


zeta_set = c(0.01,0.1,1)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  data = simu.change.regression(d0, cpt_true, p, n, sigma, kappa)
  X = data$X
  y = data$y
  cpt_true = data$cpt_true
  result = CV.search.DP.LR.regression(y, X, gamma_dp_set, lambda_dp_set, zeta_set, delta, eps = 0.001)
  min_zeta_idx = 1 + which.min(unlist(result$test_error)) %/% (length(gamma_dp_set) * length(lambda_dp_set))
  min_idx = as.vector(arrayInd(which.min(result$test_error[[min_zeta_idx]]), dim(result$test_error[[min_zeta_idx]])))
  cpt_est_lr[[time]] = unlist(result$cpt_hat[[min_zeta_idx]][min_idx[1], min_idx[2]])
  error_lr[time] = Hausdorff.dist(cpt_est_lr[[time]], cpt_true)
  len_est_lr[time] = unlist(result$K_hat[[min_zeta_idx]][min_idx[1], min_idx[2]])
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
cpt_est_lr ## estimated change points
error_lr ##Hausdorff
len_est_lr ## number of estimated change points


