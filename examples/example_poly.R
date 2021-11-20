

###########################################################
r = 1
init_coef_vec = c(1, 3)
cpt_vec = c(100, 200)
effect_size = 100
n = 300
kappa_mat = cbind(c(1, -6), c(-1, 6))
rho_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
sigma = 1
temp_mat = rho_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_mat[l+1,])^(1/(2*l+1))
}


y = gen.piece.poly(init_coef_vec, cpt_vec, kappa_mat, n, sigma)
plot.ts(y)
gamma_set = 5:19
DP_result = CV.search.DP.poly(y, r = 1, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_init = unlist(DP_result$cpt_hat[min_idx])
cpt_init
local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)


##############################################################
# Simulation 1.1 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_111.RData")



##############################################################
# Simulation 1.2 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_110.RData")


##############################################################
# Simulation 1.3 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_101.RData")




##############################################################
# Simulation 1.4 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_100.RData")





##############################################################
# Simulation 1.5 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_011.RData")




##############################################################
# Simulation 1.6 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_010.RData")




##############################################################
# Simulation 1.7 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n300_equalspace_001.RData")




##############################################################
# Simulation 2.1 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))

##############################################################
# Simulation 2.2 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))



##############################################################
# Simulation 2.3 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))


##############################################################
# Simulation 2.4 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))

##############################################################
# Simulation 2.5 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n


sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))


##############################################################
# Simulation 2.6 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))


##############################################################
# Simulation 2.7 New
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 100)
effect_size = 50
n = 150
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
mean(error_init)/n
mean(error_lr)/n


sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))





##############################################################
# Simulation 3.1 NEW
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_111.RData")


##############################################################
# Simulation 3.2
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n


sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_110.RData")



##############################################################
# Simulation 3.3
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n


sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_101.RData")



##############################################################
# Simulation 3.4
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_011.RData")



##############################################################
# Simulation 3.5
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size

pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_100.RData")



##############################################################
# Simulation 3.6
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 9, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))


save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_010.RData")



##############################################################
# Simulation 3.7
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(0, 0, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max) - 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min) - 1
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1]+1,1], rho_kl_mat[r_k[2]+1,2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = 2, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_equalspace_001.RData")




##################################################################################
# Cubic functions
##############################################################
# Simulation 4.1
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1111.RData")


##############################################################
# Simulation 4.2
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1110.RData")



##############################################################
# Simulation 4.3
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1101.RData")


##############################################################
# Simulation 4.4
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1011.RData")



##############################################################
# Simulation 4.5
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1100.RData")



##############################################################
# Simulation 4.6
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1010.RData")



##############################################################
# Simulation 4.7
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1001.RData")



##############################################################
# Simulation 4.8
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_1000.RData")



##############################################################
# Simulation 4.9
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0111.RData")


##############################################################
# Simulation 4.10
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0110.RData")


##############################################################
# Simulation 4.11
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0101.RData")



##############################################################
# Simulation 4.12
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0100.RData")


##############################################################
# Simulation 4.13
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0011.RData")


##############################################################
# Simulation 4.14
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0010.RData")


##############################################################
# Simulation 4.15
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(100, 200)
effect_size = 100
n = 300
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n300_equalspace_0001.RData")



##############################################################
# Simulation 5.1
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1111.RData")


##############################################################
# Simulation 5.2
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1110.RData")



##############################################################
# Simulation 5.3
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1101.RData")


##############################################################
# Simulation 5.4
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1011.RData")



##############################################################
# Simulation 5.5
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, -9, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1100.RData")



##############################################################
# Simulation 5.6
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1010.RData")



##############################################################
# Simulation 5.7
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1001.RData")



##############################################################
# Simulation 5.8
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(3, 0, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_1000.RData")



##############################################################
# Simulation 5.9
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0111.RData")


##############################################################
# Simulation 5.10
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0110.RData")


##############################################################
# Simulation 5.11
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0101.RData")



##############################################################
# Simulation 5.12
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, -9, 0, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0100.RData")


##############################################################
# Simulation 5.13
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 27, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0011.RData")


##############################################################
# Simulation 5.14
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 27, 0)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0010.RData")


##############################################################
# Simulation 5.15
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 3
init_coef_vec = c(2, 2, 9, 27)
cpt_true = c(150, 300)
effect_size = 150
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(-3, -9, 27, -81), c(0, 0, 0, -81)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))


## obtain signal strength
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
rho_k = apply(rho_kl_mat, MARGIN = 2, max)
l_k = apply(rho_kl_mat, MARGIN = 2, which.max)
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_true)*effect_size^(2*r_k+1)*sigma^2*log(n)/c(rho_kl_mat[r_k[1],1], rho_kl_mat[r_k[2],2]))^(1/(2*r_k+1))/effect_size



pb = txtProgressBar(min = 0, max = RR, style = 3)
counter = 0
for (time in 1:RR){
  set.seed(123+time)
  y = gen.piece.poly(init_coef_vec, cpt_true, kappa_mat, n, sigma)
  #plot.ts(y)
  DP_result = CV.search.DP.poly(y, r = r, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_est_init[[time]] = unlist(DP_result$cpt_hat[min_idx])
  error_init[time] = Hausdorff.dist(c(0, cpt_est_init[[time]], n), c(0, cpt_true, n))
  len_est_init[time] = unlist(DP_result$K_hat[min_idx])
  if(length(cpt_est_init[[time]]) == 0){
    cpt_est_lr[[time]] = cpt_est_init[[time]]
    error_lr[time] = error_init[time]
    len_est_lr[time] = len_est_init[time]
  }else{
    cpt_est_lr[[time]] = local.refine.poly(cpt_est_init[[time]], y, r = r, delta_lr = 5)
    error_lr[time] = Hausdorff.dist(c(0, cpt_est_lr[[time]], n), c(0, cpt_true, n))
    len_est_lr[time] = length(cpt_est_lr[[time]])
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
}
#local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)
scaled_Hdist_init = mean(error_init)/n
scaled_Hdist_lr = mean(error_lr)/n

save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r3_n450_equalspace_0001.RData")








r = 2
init_coef_vec = c(0, 0, 0)
cpt_vec = c(50, 100)
effect_size = 50
n = 150
kappa_mat = cbind(c(2, -5, 5), c(2, 5, -5))
rho_kl_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
sigma = 1
## obtain r_k
temp_mat = rho_kl_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_kl_mat[l+1,])^(1/(2*l+1))
}
r_k_vec = apply(temp_mat, MARGIN = 2, which.min)
(length(cpt_vec)*effect_size^(2*r_k_vec+1)*sigma^2*log(n)/c(rho_kl_mat[r_k_vec[1],1], rho_kl_mat[r_k_vec[2],2]))^(1/(2*r_k_vec+1))

y = gen.piece.poly(init_coef_vec, cpt_vec, kappa_mat, n, sigma)
plot.ts(y)
gamma_set = c(0.01, 0.1, 0.5, 1, 3:19)
DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_init = unlist(DP_result$cpt_hat[min_idx])
cpt_init
local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)




r = 3
init_coef_vec = c(1, 3, -5, 10)
cpt_vec = c(100, 200)
effect_size = 100
n = 300
kappa_mat = cbind(c(2, -10, 20, 30), c(-2, 10, -20, -30))
rho_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
sigma = 1
temp_mat = rho_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_mat[l+1,])^(1/(2*l+1))
}
  

y = gen.piece.poly(init_coef_vec, cpt_vec, kappa_mat, n, sigma)
plot.ts(y)
gamma_set = 5:19
DP_result = CV.search.DP.poly(y, r = 3, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_init = unlist(DP_result$cpt_hat[min_idx])
cpt_init
local.refine.poly(cpt_init, y, r = 3, delta_lr = 5)

