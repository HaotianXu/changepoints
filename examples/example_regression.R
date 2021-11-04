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
len_est = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error = rep(0, RR)
nrowmat = length(gamma_dp_set)
ncolmat = length(lambda_dp_set)
for (time in 1:RR){
  data = simu.change.regression(d0, change.point, p, n, sigma, kappa)
  X = data$X
  y = data$y
  cpt_true = data$cpt_true
  result = CV.search.DP.regression(y, X, gamma_dp_set, lambda_dp_set, delta, eps = 0.001)
  ##output a table which contains estimated change points, number of estimated change points, 
  ##validation loss and training error.
  min_idx = as.vector(arrayInd(which.min(result$test_error), dim(result$test_error)))
  cpt_est_init[[time]] = unlist(result$cpt_hat[min_idx[1], min_idx[2]])
  error[time] = Hausdorff.dist(cpt_est_init[[time]], cpt_true)
  len_est[time] = unlist(result$K_hat[min_idx[1], min_idx[2]])
  
  
}
cpt_est_init ## estimated change points
error ##Hausdorff
len_est ## number of estimated change points




