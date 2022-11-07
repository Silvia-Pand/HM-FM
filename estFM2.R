estFM2 <- function(Y,l,Mu0=NULL,Si0=NULL,la0=NULL,Rho0=NULL,VV0=NULL,Z=NULL,piv0 = NULL,
                   rand.ini=FALSE,nrep=0){

# fit constrained finite mixture model with constraints corresponding to an HM model for multivariate
# normal data
#
# INPUT:
# Y        = array of data of dimension n x r x T
# l        = number of latent states
# Mu0      = initial value of Mu
# Si0      = initial value of Si
# la0      = initial value of la
# Rho0     = initial value of Rho 
# VV0      = matrix of latent configurations
# Z        = design matrix for piv
# piv0     = initial value of piv
# rand.ini = for random initial assignments of individuals to latent configurations
# nrep     = if greater than 0 the estimation is repeated nrep*l times from random starting
#
# OUTPUT:
# lk = final log-likelihood value
# Mu = estimate of mean parameters
# Si = estimate of variance-covariance matrix
# Lc = global decoding

# ---- preliminaries ----
  n = dim(Y)[1]
  r = dim(Y)[2]
  TT = dim(Y)[3]
  
# ---- if there is only one state ----
  if(l==1){
    Y1 = matrix(aperm(Y,c(1,3,2)),n*TT,r)
    Mu = t(colMeans(Y1))
    Tmp = Y1-rep(1,n*TT)%*%Mu
    Si =  (t(Tmp)%*%Tmp)/(n*TT)
    lk = sum(dmvnorm(Y1,Mu,Si,log=TRUE))
    np = r+r*(r+1)/2
    aic = -2*lk+2*np
    bic = -2*lk+log(n)*np
    out = list(lk=lk,Mu=Mu,Si=Si,aic=aic,bic=bic)
    return(out)
  }

# ---- matrix of latent configurations ----
  if(is.null(VV0)) VV = sqg(rep(l-1,TT))+1 else VV = VV0
  k = nrow(VV)

# ---- starting values ----
  if(rand.ini){
    # Pp = matrix(0,n,k)
    # Pp[cbind(1:n,sample(1:k,n,rep=TRUE))] = 1
    # Pp1 = array(0,c(n,TT,l))
    # if(is.null(VV0)){
    #   for(t in 1:TT) for(v in 1:l) Pp1[,t,v] = rowSums(Pp[,VV[,t]==v])
    # }else{
    #   for(t in 1:TT) for(v in 1:l){
    #     tmp = VV[,t]==v
    #     if(sum(tmp)>0) Pp1[,t,v] = rowSums(Pp[,tmp,drop=FALSE])
    #   }
    # }
    # tot = apply(Pp1,3,sum)
# M-step
    # piv = colSums(Pp)/n
    # if(!is.null(Z)){
    #   be = c(t(Z)%*%piv)/colSums(Z)
    #   piv = c(Z%*%be)
    # }
    # Mu = matrix(0,l,r)
    # for(t in 1:TT) Mu = Mu+t(Pp1[,t,])%*%Y[,,t]
    # Mu = (1/tot)*Mu
    # Si = matrix(0,r,r)
    # for(t in 1:TT) for(v in 1:l){
    #   Tmp = Y[,,t]-rep(1,n)%o%Mu[v,]
    #   Si= Si+t(Tmp)%*%(Pp1[,t,v]*Tmp)/(n*TT)
    # }
    Y1 = matrix(aperm(Y,c(1,3,2)),n*TT,r)
    piv = runif(k); piv = piv/sum(piv)
    Mu = rmvnorm(k,colMeans(Y1),cov(Y1))
    Si = cov(Y1)*2*runif(1)
  }else{
    if(is.null(piv0)){
      if(is.null(la0)){
        piv = rep(1/k,k)
      }else{
        piv = rep(1,k)
        for(j in 1:k){
          piv[j] = la0[VV[j,1]]
          for(t in 2:TT) piv[j] = piv[j]*Rho0[VV[j,t-1],VV[j,t],t]
        }
        piv = piv/sum(piv)
      }
    }else{
      piv = piv0
    }
    if(is.null(Mu0) || is.null(Si0)) Y1 = matrix(aperm(Y,c(1,3,2)),n*TT,r)
    if(is.null(Mu0)) Mu = apply(Y1,2,quantile,probs=(1:k)/(k+1)) else Mu=Mu0
    if(is.null(Si0)) Si = cov(Y1) else Si=Si0
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
    if(is.null(VV0)){
      for(t in 1:TT) for(v in 1:l) Pp1[,t,v] = rowSums(Pp[,VV[,t]==v])
    }else{
      for(t in 1:TT) for(v in 1:l){
        tmp = VV[,t]==v
        if(sum(tmp)>0) Pp1[,t,v] = rowSums(Pp[,tmp,drop=FALSE])
      }
    }
    tot = apply(Pp1,3,sum)
# M-step
    piv = colSums(Pp)/n
    if(!is.null(Z)){
      be = c(t(Z)%*%piv)/colSums(Z)
      piv = c(Z%*%be)
    }
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
  dimnames(Mu) = list(state=1:l,variable=1:r)

# ---- global decoding ----
  lPjm = apply(lPj,1,max)
  lPj1 = lPj - lPjm
  Pp = exp(lPj1); Pp = (1/rowSums(Pp))*Pp
  Lc1 = apply(Pp,1,which.max)
  Lc = matrix(0,n,TT)
  for(i in 1:n) Lc[i,] = VV[Lc1[i],]

# ---- model selection ----
  if(is.null(Z)){
    np = r*l+r*(r+1)/2+k-1
  }else{
    np = r*l+r*(r+1)/2+ncol(Z)-1
  }
  aic = -2*lk+2*np
  bic = -2*lk+log(n)*np
  out = list(lk=lk,Mu=Mu,Si=Si,piv=piv,Lc=Lc,np=np,aic=aic,bic=bic,VV=VV)

# ---- repeat estimation from random starting values if necessary ----
  if(nrep>0){
    lkv = rep(0,nrep*l+1)
    lkv[1] = lk
    for(it in 1:(nrep*l)){
      cat("\n")
      cat("repetition from random statrting values =",it,"/",nrep*l,"\n")
      outr <- estFM2(Y,l,rand.ini=TRUE)
      lkv[it+1] = outr$lk
      if(outr$lk>out$lk) out = outr
    }
    out$lkv = lkv
  }

# ---- output ----
  return(out)

}