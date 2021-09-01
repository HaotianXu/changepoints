d0 = 5# the number of nonzero elements = p-sparsity
RR = 1## replication
sigma = 1## variance of data
kappa = 5## spectral norm
delta = 5 ## minimal spacing to help DP more robust
n = 200## sample size, need to be the multiple of 2*gridsize
p = 50 ## dimensionality
change.point = c(80,170)
gamma.dp.set = c(0.01,0.5,1,5,10,50)
lambda.dp.set = c(0.01,0.1,1,1.5)#c(0.01,0.1,1)#c(0.01,0.1,1)
len.est = matrix(0,RR,1)## number of estimated change points
init.change.point = matrix(0,20,RR) ## estimated change points
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
for (time in 1:RR){
  data = simu.change.regression(d0, change.point, p, n, sigma, kappa)
  X = data$X
  y = data$y
  true = data$cpt.true
  result = CV.search.DP.regression(y, X, gamma.dp.set, lambda.dp.set, delta)
  ##output a table which contains estimated change points, number of estimated change points, 
  ##validation loss and training error.
  min_idx = as.vector(arrayInd(which.min(result$test_error), dim(result$test_error)))
  selected.change.point = result$cpt_hat[min_idx[1], min_idx[2]]
  #selected.change.point = part2local(DP.regression(y, X, gamma.dp.set[min_idx[2]], lambda.dp.set[min_idx[1]], delta)$partition)
  init.change.point[1:length(selected.change.point), time] = selected.change.point
  #write.csv(init.change.point, file = "kappa_5_dp_d_0_15.csv")
  error[time,1] = Hausdorff.dist(selected.change.point,change.point)
  #write.csv(error, file = "kappa_5_dp_error_d_0_15.csv")
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]
  #write.csv(len.est, file = "kappa_5_dp_ncp_d_0_15.csv")
}
init.change.point ## estimated change points
error ##Hausdorff
len.est## number of estimated change points
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])



