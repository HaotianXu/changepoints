d0 = 15# the number of nonzero elements = p-sparsity
RR = 1## replication
grid = 1##gridsize to save computational cost for DP but might increase the error, default = 1
sigma = 1## variance of data
kappa = 5## spectral norm
delta = 5 ## minimal spacing to help DP more robust
n = 200## sample size, need to be the multiple of 2*gridsize
p = 50 ## dimensionality
change.point = 100
gamma.dp.set = c(0.5,1,3,20,50,70)#c(0.001,0.005,0.01)*((ncp+1)* (sigma^2) * (d0^2))
lambda.dp.set = c(0.0001,0.001,0.01,0.1,0.3,0.8)#c(0.05,0.1,0.5)*sigma*sqrt(d0)
data =simu.change.regression(d0, change.point, p,n,sigma,kappa)
X = data$X
y = data$y
true = data$cpt.true

CV.search.DP.regression(lambda.dp.set, gamma.dp.set, y, X, delta)
