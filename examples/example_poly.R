coef.transformed = function(coef_vec, cpt1, cpt2, n){
  if(length(coef_vec) == 2){
    A = matrix(c(1, (cpt2-cpt1)/n, 0, 1), nrow = 2, byrow = TRUE)
  }else if(length(coef_vec) == 3){
    A = matrix(c(1, (cpt2-cpt1)/n, ((cpt2-cpt1)/n)^2, 0, 1, 2*(cpt2-cpt1)/n, 0, 0, 1), nrow = 3, byrow = TRUE)
  }else if(length(coef_vec) == 4){
    A = matrix(c(1, (cpt2-cpt1)/n, ((cpt2-cpt1)/n)^2, ((cpt2-cpt1)/n)^3, 0, 1, 2*(cpt2-cpt1)/n, 3*((cpt2-cpt1)/n)^2, 0, 0, 1, 3*((cpt2-cpt1)/n), 0, 0, 0, 1), nrow = 4, byrow = TRUE)
  }else{
    stop("Currently, only the linear, quadratic functions and cubic functions are considered.")
  }
  return(A %*% coef_vec)
}


gen.piece.poly = function(init_coef_vec, cpt_vec, kappa_mat, n, sigma = 1){
  r = length(init_coef_vec) - 1
  cpt_ext = c(cpt_vec, n)
  if(any(dim(kappa_mat) != c(length(init_coef_vec), length(cpt_vec)))){
    stop("kappa_mat is not correct")
  }
  result = sapply((1:cpt_vec[1])/n, function(x){sum((x - cpt_vec[1]/n)^(0:r) * init_coef_vec)})
  coef_vec = init_coef_vec + kappa_mat[,1]
  for(i in 1:length(cpt_vec)){
    result = c(result, sapply(((cpt_ext[i]+1):cpt_ext[i+1])/n, function(x){sum((x - cpt_vec[i]/n)^(0:r) * coef_vec)}))
    if(i <= length(cpt_vec)-1){
      coef_vec = coef.transformed(coef_vec, cpt_vec[i], cpt_vec[i+1], n) + kappa_mat[,i+1]
    }
  }
  result = result + rnorm(n, mean = 0, sd = sigma)
  return(result)
}



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




r = 2
init_coef_vec = c(1, 3, -5)
cpt_vec = c(100, 200)
effect_size = 100
n = 300
kappa_mat = cbind(c(2, -10, 20), c(-2, 10, -20))
rho_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
sigma = 1
temp_mat = rho_mat
for(l in 0:r){
  temp_mat[l,] = (sigma*log(n)/rho_mat[l,])^(1/(2*l+1))
}


y = gen.piece.poly(init_coef_vec, cpt_vec, kappa_mat, n, sigma)
plot.ts(y)
gamma_set = 5:19
DP_result = CV.search.DP.poly(y, r = 2, gamma_set, delta = 5)
min_idx = which.min(DP_result$test_error)
cpt_init = unlist(DP_result$cpt_hat[min_idx])
cpt_init
local.refine.poly(cpt_init, y, r = 2, delta_lr = 5)


r = 2
init_coef_vec = c(1, 3, -5)
cpt_vec = c(50, 100)
effect_size = 50
n = 150
kappa_mat = cbind(c(2, -10, 20), c(-2, 10, -20))
rho_mat = apply(kappa_mat^2, MARGIN = 2, function(x){x * effect_size^(2*(0:r)+1) / n^(2*(0:r))})
sigma = 1
temp_mat = rho_mat
for(l in 0:r){
  temp_mat[l+1,] = (sigma*log(n)/rho_mat[l+1,])^(1/(2*l+1))
}


y = gen.piece.poly(init_coef_vec, cpt_vec, kappa_mat, n, sigma)
plot.ts(y)
gamma_set = 5:19
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

