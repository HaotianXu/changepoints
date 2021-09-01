set.seed(0)
p = 60 #dimension of network
B = 3 #number of blocks
K = 3 #number of change point + 1
rho = 0.05 # sparsity of the network
RR = 10 #number of repetition

#error matrices
dist_mat_com = list()
#candidate for length of time series
N_candidate = c(80,100,120,140,160)*K
#start simulation
for(cc in 1:length(N_candidate)){
  dist_mat_com[[cc]] = matrix(0, nrow = 2, ncol = RR)
  #N is the length of time series
  N = N_candidate[cc]
  # spacing parameter, only used in computing errors
  n = N/K
  for(rr in 1:RR){
    print(rr)
    can_vec = sample(1:p, replace = F) #candidate vector for SBM
    con_mat_list = list() # connectivity matrix are stored in the list
    data = matrix(nrow = p*(p-1)/2, ncol = 0)
    # generate data at each of the intervals
    # The data generated are lower diagonal as the matrices are symmetry and 0 on the diagonal
    # need some control of sbm.mean i.e. set it to be deterministic
    A1 = matrix(c(0.2,1,0.2,1,0.2,0.2,0.2,0.2,0.2), nrow = 3)
    A2 = matrix(c(0.2,0.2,1,0.2,0.2,0.2,1,0.2,0.2), nrow = 3)
    for(i in 1:K){
      tem_mat = A1
      if(i %% 2 == 0){
        tem_mat = A2
      }
      con_mat_list[[i]] = rho * tem_mat # connectivity matrix
      data = cbind(data, simu.SBM(con_mat_list[[i]] , can_vec, n))
    }
    data_mat = data
    #data splitting
    data1_mat = data_mat[,seq(2,N,2)]
    data2_mat = data_mat[,seq(1,N-1,2)]
    
    #rho.hat is an estimation of sparsity
    rho_hat = 4 * mean(rowMeans(data_mat))
    population_change = (1:(K-1))*n
    #start with NBS
    temp = BS.network(data1_mat, data2_mat, s, ncol(data1_mat), delta)
    tau = p*rho_hat
    nbs = sort(2*as.vector(t(threshold.BS(temp, tau)$change_points["location"])))
    print("nbs = ")  ;print(nbs)
    #end of NBS
    
    #compute errors for NBS
    dist_mat_com[[cc]][1,rr] = ifelse(length(nbs) == 0, n,
                                      Hausdorff.dist(nbs, population_change))
    #start with LR
    #
    lr = local.refine.network(nbs, data, tau1 = sqrt(B*p*rho_hat), tau2 = Inf, w = 1/3)
    #end of LR
    #compute errors for LR
    dist_mat_com[[cc]][2,rr] = ifelse(length(lr) == 0, n,
                                      Hausdorff.dist(lr, population_change))
    print("lr = ");print(lr)
    print(rowMeans(dist_mat_com[[cc]])/(K*n) *RR/rr)
  }
}

