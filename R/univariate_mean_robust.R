#' @title Robust test for Normal mean model
#' @noRd
rtest = function(dat1, th1, th2){
  a = 1/length(dat1)*sum(stats::dnorm(dat1, mean = th1) > stats::dnorm(dat1, mean = th2))
  b = stats::pnorm((th1+th2)/2, mean = th1)
  c = stats::pnorm((th1+th2)/2, mean = th2)
  if(abs(a-b)>abs(a-c)){
    return(1)
  }else{
    return(0)
  }
}



#' @title Robust univariate mean estimation.
#' @description Compute robust univariate mean estimation.
#' @param y       A \code{numeric} vector of observations.
#' @param k       An \code{integer} representing the lag of differences for a sorted vector.
#' @param random  A \code{logical} scalar indicating if random sampling is used for sample splitting. Otherwise, sample is split according to odd and even indices.
#' @return  A \code{numeric} scalar of robust mean estimation.
#' @noRd
#' @author Mengchu Li
#' @references Prasad, Balakrishnan and Ravikumar (2019) <arXiv:1907.00927>.
rume = function(y, k, random = FALSE){
  n = length(y)
  if(random == TRUE){
    x = sort(sample(y, floor(n/2)))
  }else{
    x = sort(y[seq(2, n, 2)])
  }
  sec = rep(NA, length(x)-k)
  for (i in 1:(length(x)-k)){
    sec[i] = x[i+k] - x[i]
  }
  e = y[is.na(pmatch(y, x))]
  ans = e[(x[which.min(sec)] <= e) & (e <= x[which.min(sec)+k])]
  if(length(ans) == 0){
    return(sum(e < x[which.min(sec)])/length(e)*x[which.min(sec)] + sum(e > x[which.min(sec)+k])/length(e)*x[which.min(sec)+k])
  }else{
    return(mean(ans))
  }
}


#' @title Adversarially robust univariate mean change point detection.
#' @description Perform the adversarially robust change point detection method.
#' @param y          A \code{numeric} vector of observations.
#' @param h          An \code{integer} scalar representing block size.
#' @param block_num  An \code{integer} scalar representing number of blocks used when searching for local maximiser.
#' @param epsilon    A \code{numeric} scalar in (0,1) representing contamination proportion.
#' @param gaussian   A \code{logical} scalar representing whether to treat the inlier distribution (F) as Gaussian distribution.
#' @return  An \code{numeric} vector of estimated change point locations
#' @export
#' @author Mengchu Li
#' @references Li and Yu (2021) <arXiv:2105.10417>.
#' @examples
#' ### simulate data with contamination
#' obs_num = 1000
#' D = 2
#' noise = 0.1 # proportion of contamination
#' mu0 = 0
#' mu1 = 1
#' sd =1
#' idmixture = rbinom(obs_num/D, 1, 1-noise)
#' dat = NULL
#' for (j in 1:D){
#'   for (i in 1:(obs_num/(2*D))){
#'     if (idmixture[i] == 1){
#'       dat = c(dat,rnorm(1,mu0,sd))
#'     }
#'     else{
#'       dat = c(dat,rnorm(1,mu1/(2*noise),0))
#'     }
#'   }
#'   for (i in (obs_num/(2*D)+1):(obs_num/D)){
#'     if (idmixture[i] == 1){
#'       dat = c(dat,rnorm(1,mu1,sd))
#'     }
#'     else{
#'       dat = c(dat,rnorm(1,mu1/(2*noise)-(1-noise)*mu1/noise,0))
#'     }
#'   }
#' }
#' plot(dat)
#' ### perform ARC
#' ARC(dat,h = 120, epsilon = 0.1)
ARC = function(y, h, block_num = 1, epsilon, gaussian = TRUE){
  sp = (1-2*epsilon-0.1)/2
  s_h = h*block_num
  cusum = NULL
  for(i in (h+1):(length(y)-h)){
    cusum = c(cusum, abs(rume(y[(i-h):i], floor(sp*h)) - rume(y[i:(i+h)], floor(sp*h))))
  }
  #plot(cusum)
  localmax = NULL
  localm = NULL
  for(j in (s_h+1):(length(cusum)-s_h)){
    if(cusum[j] == max(cusum[(j-s_h):(j+s_h)])){
      localm = c(localm,cusum[j])
      localmax = c(localmax,j)
    }
  }
  if(gaussian){
    lambda = max(8*epsilon,0.6*sqrt(40*log(length(y))/h))
  }else{
    lambda = max(8*sqrt(epsilon),0.6*sqrt(40*log(length(y))/h))
  }
  candi = which(localm > lambda)
  return(localmax[candi]+h)
}



#' @title Automatic adversarially robust univariate mean change point detection.
#' @description Perform the adversarially robust change point detection method with automatic selection of the contamination proportion epsilon when treating the inliner distributions as Gaussian.
#' @param y           A \code{numeric} vector of observations.
#' @param t_dat       A \code{numeric} vector of observations that is used to select epsilon in the Huber contamination model.
#' @param guess_true  A \code{numeric} scalar representing a guess of epsilon value.
#' @param h           An \code{integer} scalar representing block size.
#' @param block_num  An \code{integer} scalar representing number of blocks used when searching for local maximiser.
#' @return  An \code{numeric} vector of estimated change point locations.
#' @export
#' @author Mengchu Li
#' @references Li and Yu (2021) <arXiv:2105.10417>.
#' @examples 
#' #' ### simulate data with contamination
#' obs_num = 1000
#' D = 2
#' noise = 0.1 # proportion of contamination
#' mu0 = 0
#' mu1 = 1
#' sd =1
#' idmixture = rbinom(obs_num/D, 1, 1-noise)
#' dat = NULL
#' for (j in 1:D){
#'   for (i in 1:(obs_num/(2*D))){
#'     if (idmixture[i] == 1){
#'       dat = c(dat,rnorm(1,mu0,sd))
#'     }
#'     else{
#'       dat = c(dat,rnorm(1,mu1/(2*noise),0))
#'     }
#'   }
#'   for (i in (obs_num/(2*D)+1):(obs_num/D)){
#'     if (idmixture[i] == 1){
#'       dat = c(dat,rnorm(1,mu1,sd))
#'     }
#'     else{
#'       dat = c(dat,rnorm(1,mu1/(2*noise)-(1-noise)*mu1/noise,0))
#'     }
#'   }
#' }
#' plot(dat)
#' ### perform aARC
#' aARC(dat, dat[1:200], h = 120)
aARC = function(y, t_dat, guess_true = 0.05, h, block_num = 1){
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
  noi = ((-0.4 + c(1,-1)*sqrt(0.4^2+4*2*cc))/4)^2
  noi = noi[which.min(abs(noi-guess_true))]
  est_cpt = ARC(y, h = h, block_num = block_num, epsilon = noi)
  return(est_cpt)
}

#### Implementation of the backward detection algorithm with bootstrapped U statistics at 5% significance level
#' @noRd
JMB = function(y){
  n = length(y)
  cumker = 0
  for (i in 1:(n-1)){
    for (j in (i+1):n) {
      cumker = cumker+sign(y[i]-y[j])
      #print(cumker)
    }
  }
  T_n = sqrt(n)*(choose(n,2))^{-1}*cumker
  return(T_n)
}

#' @noRd
JMB_B = function(y,B){
  n = length(y)
  T = NULL
  for (b in 1:B){
    bootstrap_n = rnorm(n)
    boots_stats = NULL
    for (i in 1:(n-1)){
      cumker = 0
      for (j in (i+1):n) {
        cumker = cumker+sign(y[i]-y[j])
      }
      boots_stats = c(boots_stats,cumker*bootstrap_n[i])
    }
    T_n = sqrt(n)*(choose(n,2))^{-1}*sum(boots_stats)
    T = c(T,T_n)
  }
  return(T)
}

#' @title Backward detection with a robust bootstrap change point test using U-statistics for univariate mean change.
#' @description Perform the backward detection method with a robust bootstrap change point test using U-statistics on univariate data.
#' @param y       A \code{numeric} vector of observations.
#' @param M       An \code{integer} scalar representing initial block size of the backward detection algorithm.
#' @param B       An \code{integer} scalar representing the number of bootstrapped samples. 
#' @param inter   A nuisance parameter.
#' @param inter_ini   A nuisance parameter.
#' @return  An \code{numeric} vector of estimated change point locations.
#' @export
#' @author Mengchu Li
#' @references Yu and Chen (2019) <arXiv:1904.03372>.
BD_U = function(y, M, B = 100, inter = NULL, inter_ini = TRUE){
  C = c(0,cumsum(rep(M,length(y)/M)))
  intv = cbind(C[1:(length(C)-2)], C[3:length(C)])
  if (inter_ini == TRUE){
    inter = intv
  }
  I = dim(inter)[1]
  Test = NULL
  cpt = NULL
  if (is.null(dim(inter)) == 1){
    tn = abs(JMB(y[inter[1]:inter[2]]))
    s = sum(tn>abs(JMB_B(y[inter[1]:inter[2]],B)))/B
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
  }else{
    for (i in 1:I){
      tn = abs(JMB(y[inter[i,1]:inter[i,2]]))
      Test = c(Test,tn)
    }
    t0 = min(Test)
    t00 = which.min(Test)
    d = y[inter[t00,1]:inter[t00,2]]
    s = sum(t0>abs(JMB_B(d,B)))/B
    if (s >= 0.95){
      cpt = c(cpt,mean(inter[t00,]))
      inter = inter[-t00,]
      cpt = c(cpt,BD_U(y,M,B,inter,FALSE))
      return(sort(cpt))
    }
    else if (t00 > 1 & t00 <I){
      inter[t00-1,2] = inter[t00,2]
      inter[t00+1,1] = inter[t00,1]
      inter = inter[-t00,]
      #print(inter)
      BD_U(y,M,B,inter,FALSE)
    }
    else if (t00 == 1){
      inter[1,2] = inter[2,2]
      inter = inter[-2,]
      #print(inter)
      BD_U(y,M,B,inter,FALSE)
    }
    else if (t00 == I){
      inter[I,1] = inter[I-1,1]
      inter = inter[-(I-1),]
      #print(inter)
      BD_U(y,M,B,inter,FALSE)
    }
  }
}







#' @title Robust wild binary segmentation for univariate mean change points detection.
#' @description Perform a robust version of the wild binary segmentation method using Huber loss. 
#' @param y       A \code{numeric} vector of observations.
#' @param s       A \code{numeric} scalar representing the starting index of an interval.
#' @param e       A \code{numeric} scalar representing the ending index of an interval.
#' @param Alpha   An \code{integer} vector of starting indices of random intervals.
#' @param Beta    An \code{integer} vector of ending indices of random intervals.
#' @param K       A \code{numeric} scalar representing the robustification parameter in the Huber loss.
#' @param delta   A positive \code{integer} scalar of minimum spacing.
#' @param level   Should be fixed as 0.
#' @return  An object of \code{\link[base]{class}} "BS", which is a \code{list} with the following structure:
#'  \item{S}{A vector of estimated change points (sorted in strictly increasing order).}
#'  \item{Dval}{A vector of values of CUSUM statistic based on KS distance.}
#'  \item{Level}{A vector representing the levels at which each change point is detected.}
#'  \item{Parent}{A matrix with the starting indices on the first row and the ending indices on the second row.}
#' @export
#' @author Mengchu Li & Haotian Xu
#' @references Fearnhead & Rigaill (2019) <doi:10.1080/01621459.2017.1385466>.
#' @seealso \code{\link{thresholdBS}} for obtaining change points estimation.
WBS.uni.rob = function(y, s, e, Alpha, Beta, K = 1.345, delta, level = 0){
  obs_num = length(y)
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(M == 0){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, M)
    b = rep(0, M)
    for(m in 1:M){
      r1 = y[(Alpha_new[m]+1):Beta_new[m]] - huber_mean(y[(Alpha_new[m]+1):Beta_new[m]], K)
      r1 = sapply(r1, function(x){if(abs(x) < K){return(x)}else if(x < -K){return(-K)}else{return(K)}})
      r2 = cumsum(r1)
      r3 = rep(0, length(r2)-1)
      for(i in 1:(length(r2)-1)){
        r3[i] = r2[i]^2*length(r2)/((length(r2)-i)*i)
      }
      best_value = max(r3)
      best_t = which.max(r3) + Alpha_new[m]
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = WBS.uni.rob(y, s, b[m_star]-1, Alpha, Beta, K, delta, level)
  temp2 = WBS.uni.rob(y, b[m_star], e, Alpha, Beta, K, delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}
