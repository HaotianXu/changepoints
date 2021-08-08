p = 20
sigma = 1
n = 30
v1 = 2*(seq(1,p,1)%%2) - 1
v2 = -v1
AA = matrix(0, nrow = p, ncol = p-2)
A1=cbind(v1,v2,AA)
A2=cbind(v2,v1,AA)
A3=A1
data = simu.change.VAR1(sigma, p, 2*n+1, A1)
data = cbind(data, simu.change.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
data = cbind(data, simu.change.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
N = ncol(data)
X_curr = data[,1:(N-1)]
X_futu = data[,2:N]
parti = DP.VAR1(gamma = 1, delta = 5, X_futu, X_curr, lambda = 1)$partition
localization = part2local(parti)


