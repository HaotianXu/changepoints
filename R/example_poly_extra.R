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




##############################################################
# Simulation 2.5 Extra
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
# kappa_mat = cbind(c(3, 9, -27), c(0, 404.1404, 0)) # jump sizes of coefficients for reparametrized model
# (length(cpt_true)*sigma^2*log(n))^(-1) * 3^3 * n
kappa_mat = cbind(c(3, 9, -27), c(0, 0, -54443.15)) # jump sizes of coefficients for reparametrized model
# (length(cpt_true)*sigma^2*log(n))^(-2) * 3^5 * n^2
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
mean(error_init)/n
mean(error_lr)/n