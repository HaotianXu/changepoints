##############################################################
# Simulation 1
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
spacing_vec = c(50, 50)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 100) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_50_350.RData")

##############################################################
# Simulation 2
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 150)
spacing_vec = c(50, 100)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 150) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_100_300.RData")



##############################################################
# Simulation 3
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 200)
spacing_vec = c(50, 150)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_150_250.RData")


##############################################################
# Simulation 4
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 250)
spacing_vec = c(50, 200)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 250) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_200_200.RData")


##############################################################
# Simulation 5
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 300)
spacing_vec = c(50, 150)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_250_150.RData")



##############################################################
# Simulation 6
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 350)
spacing_vec = c(50, 100)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 350) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_300_100.RData")

##############################################################
# Simulation 7
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(50, 400)
spacing_vec = c(50, 50)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 50) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 400) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_50_350_50.RData")



##############################################################
# Simulation 8
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
spacing_vec = c(100, 100)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_100_100_250.RData")


##############################################################
# Simulation 9
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 250)
spacing_vec = c(100, 150)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 250) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_100_150_200.RData")


##############################################################
# Simulation 10
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 300)
spacing_vec = c(100, 150)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 300) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_100_200_150.RData")



##############################################################
# Simulation 11
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(100, 350)
spacing_vec = c(100, 100)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 350) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_100_250_100.RData")


##############################################################
# Simulation 12
RR = 100
len_est_init = rep(0, RR)## number of estimated change points
cpt_est_init = vector("list", RR) ## estimated change points
error_init = rep(0, RR)
len_est_lr = rep(0, RR)## number of estimated change points
cpt_est_lr = vector("list", RR) ## estimated change points
error_lr = rep(0, RR)
r = 2
init_coef_vec = c(-2, 2, 9)
cpt_true = c(200, 250)
spacing_vec = c(50, 50)
effect_size = min(spacing_vec)
n = 450
sigma = 1
gamma_set = c(0.01, 0.05, 0.1, 0.5, 1:16)
kappa_mat = cbind(c(3, 9, -27), c(-3, 9, -27)) # jump sizes of coefficients for reparametrized model

init_coef_vec+kappa_mat[,1] # coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n) # reparametrized coefficients after the first changepoint
coef.repara(init_coef_vec+kappa_mat[,1], cpt_true[1], cpt_true[2], n)+kappa_mat[,2] # coefficients after the second changepoint
plot.ts(gen.piece.poly.noiseless(init_coef_vec, cpt_true, kappa_mat, n, sigma))

## obtain signal strength
rho_kl_mat = kappa_mat^2 * rbind(spacing_vec, spacing_vec^3/n^2, spacing_vec^5/n^4) 
row.names(rho_kl_mat) <- NULL
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

sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 200) <= 3)}))
sum(sapply(1:length(cpt_est_init), function(t){sum(abs(cpt_est_init[[t]] - 250) <= 3)}))



save(r, init_coef_vec, cpt_true, effect_size, n, sigma, kappa_mat, cpt_est_init, error_init, len_est_init, cpt_est_lr, error_lr, len_est_lr, scaled_Hdist_init, scaled_Hdist_lr, file = "C:/Users/haotian/Documents/GitHub/changepoints/examples/poly_r2_n450_changespace_200_50_200.RData")


