NBS = function(Y, delta, s, e, N){
  S = NULL
  Dval = NULL
  if(e-s <= delta){
    S = c(S, NULL)
    Dval = c(Dval, NULL)
    return(list(S = S, Dval = Dval))
  }else{
    a = rep(0, e-s-1)
    for(t in (s+1):(e-1)){
      a[t-s] = Delta_se_t(Y, s, e, t, N)
    }
    best_value = max(a)
    best_t = which.max(a) + s
    S = c(S, best_t)
    Dval = c(Dval, best_value)
    temp1 = NBS(Y, delta, s, best_t-1, N)
    temp2 = NBS(Y, delta, best_t, e, N)
    S = c(temp1$S, S, temp2$S)
    Dval = c(temp1$Dval, Dval, temp2$Dval)
    return(list(S = S, Dval = Dval))
  }
}


points(x = tail(temp$S[order(temp$Dval)],20), y =rep(0, 20), col = "red")


Delta_se_t = function(Y, s, e, t, N){
  T =  dim(Y)[2]
  n =  dim(Y)[1]
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te = sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux = Y[, s:t]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y = Y[, s:e]
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st = temp(vec_y)# temp(grid)
  
  aux = Y[, (t+1):e]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te = temp(vec_y)# temp(grid)
  
  temp = sqrt(n_st * n_te / n_se) * max(abs(Fhat_te - Fhat_st))
  return(temp)
}




NBS_full =  function(y, z, gam, N){
  n = max(N)
  T = dim(y)[1]
  temp1 = new_BS(z, gam,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-5})
  tau_grid =  tau_grid[which(is.na(tau_grid)==FALSE)] ### *
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    if(length(aux)==0)
    {
      break;
    }
    S[[j]] = sort(aux)
  }
  #}
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  
  lamb =log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:length(S))#)
  {
    if(length(S[[j]])==0)
    {
      j = j+1;
      
      if(j>length(S))
        break;
    }
    
    B2  =  S[[j]]
    if(j==length(S))
    {
      B1 = NULL
    }
    if(j< length(S))
    {
      B1 = S[[j+1]]
    }
    temp = setdiff(B2,B1)
    
    st =  -10^15
    #Delta_se_t(z,eta1,eta2,eta,N)^2 
    for(l in 1:length(temp))
    {
      eta =  temp[l]
      
      if(j == length(S))
      {
        eta1 = 1
        eta2 = T
      }
      if(j < length(S))
      {
        for(k in 1:length(S[[j+1]]))
        {
          if(S[[j+1]][k]> eta  )
            break;
        }
        if(S[[j+1]][k]> eta )
        {
          eta2 = S[[j+1]][k]
          
          if(k ==1)
            eta1 = 1
          
          if(k > 1)
            eta1 = S[[j+1]][k-1]+1
        }
        if(S[[j+1]][k]< eta )
        {
          eta1 = S[[j+1]][k]+1
          eta2 = T
        }
      }
      st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
      # print(st_au)
      if(st_aux> st)
      {
        st = st_aux
      }
    }###  close for defining  eta1 and eta2
    
    
    # print(c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb)
    if(st >   lamb)
    {
      #B1 = B2
      return(B2)
    }
    # print(st)
  }
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb
  
  return(B1)
  # 
  
}