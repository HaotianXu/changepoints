set.seed(123)

#START define parameters
rho_can = c(0.1, 0.4)#candidate of rho (rho can not be too large so that the time series is stable)
ll = length(rho_can)
RR = 3# 100 # number of repetition
p = 20 # dimension
sigma = 1 # standard deviation of error terms
delta1 = 10 # 10
delta.local = 10 # 10
# construct true transition matrices
n = 30
v1 = 2*(seq(1,p,1)%%2) - 1
v2 = -v1 
AA = matrix(0, nrow = p, ncol = p-2)
A1=cbind(v1,v2,AA)*0.1 # transition matrix for the first segment
A2=cbind(v2,v1,AA)*0.1 # transition matrix for the second segment
A3=A1 # transition matrix for the third segment
# generate data
data = simu.VAR1(sigma, p, 2*n+1, A1)
data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
#FINISH generate data
N = ncol(data) # sample size
X_curr = data[,1:(N-1)]
X_futu = data[,2:N]
parti = DP.VAR1(X_futu, X_curr, gamma = 1, lambda = 1, delta = delta1)$partition
localization = part2local(parti)

#records of change point estimation
dp_rec = vector("list", length = ll)
lr_rec = vector("list", length = ll)
#checking the errors of DP and LR 
haus_dp = rep(0,ll)
haus_lr = rep(0,ll)


for(candidate in 1:ll){
  print(paste("candidate=" ,candidate))
  rho = rho_can[candidate]
  A1 = A1*rho
  A2 = A2*rho
  A3 = A3*rho

  for(rr in 1:RR){
    print(paste("rr=" ,rr))
    data = simu.VAR1(sigma, p, 2*n+1, A1)
    data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
    data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
    #For convenient due to data splitting, col(data) (or the length of the time series) is odd.
    #the change points are at 2*n and 4*n
    lambda_lasso_set = c(3,4,10,20)
    gamma_set = seq(1, 50,2) 
    #START: CV for DP
    start.time <- Sys.time()
    dp_result = CV.search.DP.VAR1(data, gamma_set, lambda_lasso_set, delta1)
    end.time <- Sys.time(); time.taken <- end.time - start.time; print(time.taken)
    min_idx = as.vector(arrayInd(which.min(dp_result$test_error), dim(dp_result$test_error)))
    X_curr = data[,1:(N-1)]
    X_futu = data[,2:N]
    dp_estimate = part2local(DP.VAR1(X_futu, X_curr, gamma = gamma_set[min_idx[2]], lambda = lambda_lasso_set[min_idx[1]], delta = delta1)$partition)
    dp_estimate_ext = c(0, dp_estimate, 6*n)
    print("dp=")
    print(dp_estimate_ext)
    dp_rec[[candidate]][[rr]] = list(dp_estimate_ext)
    haus_dp[candidate] = haus_dp[candidate] + Hausdorff.dist(dp_estimate_ext, c(0, 2*n, 4*n, 6*n))
    print(haus_dp/rr)
    #END:  DP
    
    zeta_group_set = c(0.5, 1, 1.5)
    #START  local refinement#include X and Y
    lr_estimate = local.refine.CV.VAR1(dp_estimate, data, zeta_group_set, delta.local, w = 1/3)
    lr_estimate_ext = c(0, lr_estimate, 6*n)
    print( "lr=")
    print(lr_estimate_ext)
    lr_rec[[candidate]][[rr]] = list(lr_estimate)
    haus_lr[candidate] = haus_lr[candidate] + Hausdorff.dist(lr_estimate_ext, c(0, 2*n,4*n,6*n) )
    print(haus_lr/rr)
    #END local refinement
  }
}











