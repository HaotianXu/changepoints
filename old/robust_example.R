
################################ example - hide

obs_num = 1000
D = 2

noise = 0.1
mu0 = 0
mu1 = 1
sd =1


idmixture = rbinom(obs_num/D,1,1-noise)
dat = NULL

for (j in 1:D){
  for (i in 1:(obs_num/(2*D))){
    if (idmixture[i] == 1){
      dat = c(dat,rnorm(1,mu0,sd))
    }
    else{
      dat = c(dat,rnorm(1,mu1/(2*noise),0))
    }
  }
  
  for (i in (obs_num/(2*D)+1):(obs_num/D)){
    if (idmixture[i] == 1){
      dat = c(dat,rnorm(1,mu1,sd))
    }
    else{
      dat = c(dat,rnorm(1,mu1/(2*noise)-(1-noise)*mu1/noise,0))
    }
  }
}

plot(dat)
ARC(dat,h = 120,epsilon = 0.1)
aARC(dat,dat[1:200],h = 120)
BD_U(dat,100)
intervals = WBS.intervals(M = 300, lower = 1, upper = length(dat))
temp = WBS.uni.rob(dat, 1, length(dat), intervals$Alpha, intervals$Beta, delta = 5)
thresholdBS(temp, 15)
########################################### example - fake

idmixture = rbinom(obs_num/D,1,1-noise)
dat = NULL
for (j in 1:D){
  for (i in 1:(obs_num/(2*D))){
    if (idmixture[i] == 1){
      dat = c(dat,rnorm(1,mu0,sd))
    }
    else{
      dat = c(dat,rnorm(1,-3,0))
    }
  }
  
  for (i in (obs_num/(2*D)+1):(obs_num/D)){
    if (idmixture[i] == 1){
      dat = c(dat,rnorm(1,mu0,sd))
    }
    else{
      dat = c(dat,rnorm(1,3,0))
    }
  }
}
plot(dat)

ARC(dat,h = 120,epsilon = 0.1)
aARC(dat,dat[1:200],h = 120)
BD_U(dat,100)
intervals = WBS.intervals(M = 300, lower = 1, upper = length(dat))
temp = WBS.uni.rob(dat, 1, length(dat), intervals$Alpha, intervals$Beta, delta = 5)
thresholdBS(temp, 15)
