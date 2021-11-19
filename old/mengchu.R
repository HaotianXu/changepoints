library(pracma)
#library(devtools)
#install_github("guillemr/robust-fpop")
#require(robseg)
library(wbs)
library(qrmix)


############################ ARC and aARC algorithms 

#' @title ...
#' @noRd
result <- function(a, b, c){
  delta = b^2-4*a*c
  if(delta > 0){ # first case D>0
    x_1 = (-b+sqrt(delta))/(2*a)
    x_2 = (-b-sqrt(delta))/(2*a)
    result = c(x_1,x_2)
  }
  else if(delta == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

#' @title Robust test for Normal mean model
#' @noRd
rtest = function(dat1, th1, th2){
  a = 1/length(dat1)*sum(dnorm(dat1, mean = th1)>dnorm(dat1, mean = th2))
  b = pnorm((th1+th2)/2, mean = th1)
  c = pnorm((th1+th2)/2, mean = th2)
  if (abs(a-b)>abs(a-c)){
    return(1)
  }
  else{
    return(0)
  }
}



#' @title Robust Univariate Mean Estimation
#' @description Compute robust univariate mean estimation
#' @param y  A \code{numeric} vector of observations
#' @param k  An \code{integer} representing the lag of differences for a sorted vector
#' @return  A numeric scalar of robust mean estimation.
#' @export
#' @author Mengchu Li
#' @references Prasad, Adarsh, Sivaraman Balakrishnan, and Pradeep Ravikumar. "A unified approach to robust mean estimation." arXiv preprint arXiv:1907.00927 (2019).
RUME = function(y, k){
  n = length(y)
  x = sort(sample(y, floor(n/2)))
  sec = rep(NA, length(x)-k)
  for (i in 1:(length(x)-k)){
    sec[i] = x[i+k] - x[i]
  }
  e = y[is.na(pmatch(y, x))]
  ans = e[(x[which.min(sec)] <= e) & (e <= x[which.min(sec)+k])]
  return(mean(ans))
}


#' @title Adversarially robust change point detection
#' @description Perform the adversarially robust change point detection method.
#' @param y        A \code{numeric} vector of observations
#' @param h        An \code{integer} scalar representing block size
#' @param lambda   A \code{numeric} scalar representing the lag of differences for a sorted vector
#' @param epsilon  A \code{numeric} scalar representing ....
#' @param s_h      An \code{integer} scalar representing ...
#' @param heavyt   A \code{logical} scalar representing ...
#' @return  ....
#' @export
#' @author Mengchu Li
#' @references Li, M., & Yu, Y. (2021). Adversarially Robust Change Point Detection. arXiv preprint arXiv:2105.10417.
ARC = function(y, h = floor(20*log(length(y))), lambda = max(8*noi, 0.6), epsilon, s_h = floor(20*log(length(y))), heavyt = FALSE){
  sp = (1-2*epsilon-0.1)/2
  cusum = NULL
  for (i in (h+1):(length(y)-h)) {
    cusum = c(cusum, abs(rume(y[(i-h):i], floor(sp*h)) - rume(y[i:(i+h)], floor(sp*h))))
  }
  #plot(cusum)
  localmax = NULL
  localm = NULL
  for (j in (s_h+1):(length(cusum)-s_h)){
    if (cusum[j] == max(cusum[(j-s_h):(j+s_h)])){
      localm = c(localm,cusum[j])
      localmax = c(localmax,j)
    }
  }
  if (heavyt){
    lambda = max(8*sqrt(epsilon),0.6*sqrt(40*log(length(y))/h))
  }
  else{
    lambda = max(8*epsilon,0.6*sqrt(40*log(length(y))/h))
  }
  #print(lambda)
  candi = which(localm>lambda)
  return(localmax[candi]+h)
}



#' @title Automatic adversarially robust change point detection
#' @description Perform the adversarially robust change point detection method.
#' @param dat         A \code{numeric} vector of observations
#' @param t_dat       ....
#' @param guess_true  ....
#' @param h           An \code{integer} scalar representing block size
#' @param s_h         An \code{integer} scalar representing ...
#' @return  ....
#' @export
#' @author Mengchu Li
#' @references Li, M., & Yu, Y. (2021). Adversarially Robust Change Point Detection. arXiv preprint arXiv:2105.10417.
aARC = function(dat, t_dat, guess_true = 0.05, h, s_h){
  eplison = seq(0.25, 0.45, 0.001)
  est_e = NULL
  for (e in 1:length(eplison)){
    holder <- 0
    for (i in 1:5){
      holder = holder + rume(t_dat,floor(eplison[e]*length(t_dat)))
    }
    est_e = c(est_e,holder/5)
  }
  ho = rep(0, 201)
  for (e1 in 1:length(est_e)){
    for (e2 in 1:length(est_e)){
      a = min(est_e[e1],est_e[e2])
      b = max(est_e[e1],est_e[e2])
      if (a == est_e[e1]){
        ho[e1] = ho[e1] + rtest(t_dat,a,b)
      }else if (b == est_e[e1]){
        ho[e1] = ho[e1] + abs(1-rtest(t_dat,a,b))
      }
    }
  }
  ind = which.min(ho)
  eps = eplison[ind]
  cc = 0.96 - 2*eps
  noi = result(2, 0.4, -cc)^2
  noi = noi[which.min(abs(noi-guess_true))]
  est_cpt = ARC(dat,eps = noi,h = h,s_h = s_h)
  return(est_cpt)
}

#### Implementation of the backward detection algorithm with bootstrapped U statistics at 5% significance level

JMB = function(dat){
  n = length(dat)
  cumker = 0
  for (i in 1:(n-1)){
    for (j in (i+1):n) {
      cumker = cumker+sign(dat[i]-dat[j])
      #print(cumker)
    }
  }
  T_n = sqrt(n)*(choose(n,2))^{-1}*cumker
  return(T_n)
}

JMB_B = function(dat,B){
  n = length(dat)
  T = NULL
  for (b in 1:B){
    bootstrap_n = rnorm(n)
    boots_stats = NULL
    for (i in 1:(n-1)){
      cumker = 0
      for (j in (i+1):n) {
        cumker = cumker+sign(dat[i]-dat[j])
      }
      boots_stats = c(boots_stats,cumker*bootstrap_n[i])
    }
    T_n = sqrt(n)*(choose(n,2))^{-1}*sum(boots_stats)
    T = c(T,T_n)
  }
  return(T)
}


BD_U = function(dat, M, B = 100, inter = NULL, inter_ini = TRUE){
  C = c(0,cumsum(rep(M,length(dat)/M)))
  intv = cbind(C[1:(length(C)-2)], C[3:length(C)])
  if (inter_ini == TRUE){
    inter = intv
  }
  I = dim(inter)[1]
  Test = NULL
  cpt = NULL
  if (is.null(dim(inter)) == 1){
    tn = abs(JMB(dat[inter[1]:inter[2]]))
    s = sum(tn>abs(JMB_B(dat[inter[1]:inter[2]],B)))/B
    if (s>=0.95){
      cpt = c(cpt,mean(inter))
      return(cpt)
    }
    else{
      return(cpt)
    }
    return(cpt)
  }
  if (I < 1){
    return(cpt)
  }
  else{
    for (i in 1:I){
      tn = abs(JMB(dat[inter[i,1]:inter[i,2]]))
      Test = c(Test,tn)
    }
    t0 = min(Test)
    t00 = which.min(Test)
    d = dat[inter[t00,1]:inter[t00,2]]
    s = sum(t0>abs(JMB_B(d,B)))/B
    #browser()
    #print(s)
    if (s >= 0.95){
      cpt = c(cpt,mean(inter[t00,]))
      inter = inter[-t00,]
      cpt = c(cpt,BD_U(dat,M,B,inter,FALSE))
      return(sort(cpt))
    }
    else if (t00 > 1 & t00 <I){
      inter[t00-1,2] = inter[t00,2]
      inter[t00+1,1] = inter[t00,1]
      inter = inter[-t00,]
      #print(inter)
      BD_U(dat,M,B,inter,FALSE)
    }
    else if (t00 == 1){
      inter[1,2] = inter[2,2]
      inter = inter[-2,]
      #print(inter)
      BD_U(dat,M,B,inter,FALSE)
    }
    else if (t00 == I){
      inter[I,1] = inter[I-1,1]
      inter = inter[-(I-1),]
      #print(inter)
      BD_U(dat,M,B,inter,FALSE)
    }
  }
}







####################################### robust version of the wild binary segmentation algorithm 
HLoss = function(t,d){
  test = sum(Huber(d-t), c= 1.345)
  return(test)
}

wbs_r = function(dat, M, s = 1, e = length(dat), intervals = NULL, hat = 5*mad(dat)*log(length(dat)), K = 1.345){
  T = length(dat)
  interv = random.intervals(T,M)
  n = e-s
  test = NULL
  scpt = NULL
  cpt = NULL
  for (i in 1:M){
    #print(intervals[i,])
    inter = intersect(interv[i,1]:interv[i,2],s:e)
    if (length(inter) < 5){
      test = c(test,-1)
      scpt = c(scpt,-1)
    }
    else{
      interv[i,] = c(inter[1],inter[length(inter)])
      #print(sub)
      datint = dat[inter]
      #print(datint)
      #browser()
      opti = optim(0,HLoss,d = datint, method="Brent",lower = min(dat), upper = max(dat))
      theta = opti$par
      r1 = datint-theta
      L = length(r1)
      for (i in 1:L){
        if ((datint[i]-theta)>K){
          r1[i] = K
        }
        else if ((datint[i]-theta)<(-K)){
          r1[i] = -K
        }
      }
      r2 = cumsum(r1)
      #browser()
      for (i in 1:(L-1)){
        r2[i] = r2[i]^2*L/((L-i)*i)
      }
      r2 = r2[is.finite(r2)]
      Tmax = max(r2)
      tmax = which.max(r2)
      test = c(test,Tmax)
      #print(test)
      #plot(r2)
      scpt = c(scpt,tmax)
    }
  }
  #browser()
  #print(max(test))
  if (max(test) > hat){
    t1 = scpt[which.max(test)]+interv[which.max(test),1]
    #print(t1)
    #print(intervals[which.max(test),])
    cpt = c(cpt,t1)
    cpt = c(cpt,wbs_r(dat,M,s,(t1-1),intervals = interv))
    cpt = c(cpt,wbs_r(dat,M,(t1+1),e,intervals = interv))
    cpt = cpt[cpt!=(-1)]
    return(sort(cpt))
  }
  else{
    return(-1)
  }
}
