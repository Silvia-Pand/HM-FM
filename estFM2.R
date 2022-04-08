estFM2 <- function(Y,l,Mu=NULL,Si=NULL,piv0=NULL,Pi0=NULL){

# fit constrained finite mixture model with constraints corresponding to a HM model for multivariate
# normal data
  
# INPUT:
#
# Y = array of data of dimension n x r x T
#  l = number of latent states
#
# OUTPUT:
#
# lk = final log-likelihood value
# Mu = estimate of mean parameters
# Si = estimate of variance-covariance matrix
# Lc = global decoding

# ---- preliminaries ----
  n = dim(Y)[1]
  r = dim(Y)[2]
  TT = dim(Y)[3]
  k = l^TT
  VV = sqg(rep(l-1,TT))+1

# ---- starting values ----
  if(is.null(piv0)) piv = rep(1/k,k)
  else{
    piv = rep(1,k)
    for(j in 1:k){
      piv[j] = piv0[VV[j,1]]
      for(t in 2:TT) piv[j] = piv[j]*Pi0[VV[j,t-1],VV[j,t],t]
    }   
  }
  Y1 = matrix(aperm(Y,c(1,3,2)),n*TT,r)
  if(is.null(Mu)) Mu = apply(Y1,2,quantile,probs=(1:k)/(k+1))
  if(is.null(Si)) Si = cov(Y1)

# ---- compute log-likelihood ----
  lGG = array(0,c(n,l,TT))
  for(t in 1:TT) for(v in 1:l) lGG[,v,t] = dmvnorm(Y[,,t],Mu[v,],Si,log=TRUE)
  lFF = matrix(0,n,k)
  for(t in 1:TT) lFF = lFF+lGG[,VV[,t],t]
  lFFm = apply(lFF,1,max)
  lPj = lFF + rep(1,n)%o%log(piv)
  lFF1 = lFF - lFFm
  lpm = log(exp(lFF1)%*%piv)+lFFm
  lk = sum(lpm)
  print(lk)

# ---- iterate until convergence ----
  it = 0; lko = lk
  while((lk-lko)/abs(lko)>10^-10 || it==0){
    lko = lk; it = it+1
# E-step
    lPjm = apply(lPj,1,max)
    lPj1 = lPj - lPjm
    Pp = exp(lPj1); Pp = (1/rowSums(Pp))*Pp
    Pp1 = array(0,c(n,TT,l))
    for(t in 1:TT) for(v in 1:l) Pp1[,t,v] = rowSums(Pp[,VV[,t]==v])
    tot = apply(Pp1,3,sum)
# M-step
    piv = colSums(Pp)/n
    Mu = matrix(0,l,r)
    for(t in 1:TT) Mu = Mu+t(Pp1[,t,])%*%Y[,,t]
    Mu = (1/tot)*Mu
    Si = matrix(0,r,r)
    for(t in 1:TT) for(v in 1:l){
      Tmp = Y[,,t]-rep(1,n)%o%Mu[v,]
      Si= Si+t(Tmp)%*%(Pp1[,t,v]*Tmp)/(n*TT)
    }
    # ---- compute log-likelihood ----
    lGG = array(0,c(n,l,TT))
    for(t in 1:TT) for(v in 1:l) lGG[,v,t] = dmvnorm(Y[,,t],Mu[v,],Si,log=TRUE)
    lFF = matrix(0,n,k)
    for(t in 1:TT) lFF = lFF+lGG[,VV[,t],t]
    lFFm = apply(lFF,1,max)
    lPj = lFF + rep(1,n)%o%log(piv)
    lFF1 = lFF - lFFm
    lpm = log(exp(lFF1)%*%piv)+lFFm
    lk = sum(lpm)
    if(it%%100==0) print(c(it,lk,lk-lko))
  }
  if(it%%100>0) print(c(it,lk,lk-lko))

# ---- global decoding ----
  lPjm = apply(lPj,1,max)
  lPj1 = lPj - lPjm
  Pp = exp(lPj1); Pp = (1/rowSums(Pp))*Pp
  Lc1 = apply(Pp,1,which.max)
  Lc = matrix(0,n,TT)
  for(i in 1:n) Lc[i,] = VV[Lc1[i],]

# ---- output ----
  out = list(lk=lk,Mu=Mu,Si=Si,piv=piv,Lc=Lc)

}