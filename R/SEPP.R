#generate stable SEPP data

gen.seep.data= function(intercept,M,A,threshold,n,vzero){
  
  #X is the data matrix with horizontal axis being time
  if(length(vzero)==0){
    vthreshold=rpois(M, intercept)}else{
      vthreshold=vzero
      vthreshold[which(vzero>threshold)]=threshold
    }
  X=matrix(0,ncol=n,nrow=M)
  
  X[,1]= rpois(M, lambda=exp(intercept+A%*%as.matrix(vthreshold)))
  
  for ( t in 2:n){
    X.temp = X[,t-1]
    X.temp[which(X[,t-1]>threshold)]=threshold
    
    X[,t]=rpois(M,lambda= exp(intercept+A%*%X.temp))
    
    
  }
  return(X)}