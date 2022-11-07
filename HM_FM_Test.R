HM_FM_Test <- function(data=data,id = "id",time="time",lv=2:4,responses=NULL,
                       type=c("single","triplets"),al=0.05){

#Function that performs the LR test for HM-FM models
# 
# INPUT:
#data       = a data.frame in long format
#id         = a string indicating the name of the unit identifier
#time       = a string indicating the name of the time occasions
#lv         = an integer vector specifying the number of latent states for which the test has to be performed
#responses  = a character vector indicating the name of the response variables
#             If NULL colnames(data) from the third to the last column is taken
#type       = the type of test that has to be performed ("single" = the standard LR test,
#             "triplets" = the multiple test based on consecutive triplets of observations
#               based on the Bonferroni and on the Simes corrections)
#al         = significance level
#
# OUTPUT:
#lkHM    = value of the  log-likelihood of the HM model for each l in lk
#aicHM   = value of AIC index for the HM model for each l in lv
#bicHM   = value of BIC index for the HM model for each l in lv
#n       = sample size
#TT      = number of time occasions
#pvalue  = pvalue of the standard LR test for each l in lv (if type="single")   
#pvalueB  = overall pvalue of the multiple LR test for each l in lv based on the Bonferroni correction (if type="triplets")
#pvalueS  = overall pvalue of the multiple LR test for each l in lv based on the Simes correction (if type="triplets")
#Pvc     = matrix of single p-values of triples (if type="triplets")
#reject  = conclusion of the test based on al for the standard LR test (if type="single") 
#rejectB =  conclusion of the test based on al with Bonferroni rule (if type="triplets")
#rejectS =  conclusion of the test based on al with Simes rule (if type="triplets")
#selstates = selected number of latent states among those in lv
  
# ---- preliminaries ----
  type <- match.arg(type, choices = eval(formals(HM_FM_Test)$type))
  if(is.null(responses)) responses = colnames(data[,3:dim(data)[2]])
  responsesFormula <- lmestFormula(data=data,response=responses)$responsesFormula
  aicHM <-bicHM <- lkHM <- Dev1<- pvalue <- pvalueB <- pvalueS <- reject <- rejectB <- rejectS <- df <- dft<- rep(NA,max(lv))
  out <-long2matrices(id=data[,names(data)==id],time =data[,names(data)==time],Y=data[,names(data)%in%responses])
  n <- dim(out$YY)[1]
  TT <- dim(out$YY)[2]
  YY <- out$YY
  if(TT>3 & type=="triplets"){     
    Pvc = matrix(0,TT-2,max(lv))   
    dimnames(Pvc) = list(t=1:(TT-2),l=1:max(lv))   
  }  #

# ---- iterate for each l ----   
  for(l in lv){
    if(type=="single" & l^TT>200){
      cat("|--------------------------WARNING---------------------------------------|\n")
      cat("|The number of possibile sequences of latent states is particularly large|\n")
      cat("|                   Please use option triplets                           |\n")
      cat("|------------------------------------------------------------------------|\n")
    }
    cat("Estimate the HM model for l =",l,"latent states","\n")
    estHM <- lmestCont(responsesFormula = responsesFormula,
                           index = c(id,time),
                           data = data,
                           k = l,
                           modBasic = 1,
                           tol = 10^-10)
    lkHM[l] <- estHM$lk
    aicHM[l] <- estHM$aic
    bicHM[l] <- estHM$bic
    if(l>1){
      if(type=="single" || TT<=3){   
        print("Perform the standard LR test")
        cat("Estimate the FM2 model for k =",l^TT,"components","\n")
        estFM <- estFM2(aperm(YY,c(1,3,2)),l=l, Mu=t(estHM$Mu),Si=estHM$Si,la0=estHM$piv,Rho0=estHM$Pi)
        lkFM <- estFM$lk
        Dev1[l] <- -2*(lkHM[l] - lkFM)
        df[l] = l^TT-l^2
        pvalue[l] <- 1-pchisq(Dev1[l],l^TT-l^2)
        reject[l] <- pvalue[l]<al
      }else if(type=="triplets"){
        print("perform the multiple LR test")
        if(TT>3){
          lkHMc <- lkFMc <- devc <- rep(0,TT-2)
          for(t in 1:(TT-2)){
            Ytime <- data[,names(data)==time]
            Yt = data[Ytime%in%(t:(t+2)),]
            Yt$time = Yt$time-(t-1)
            cat("Estimate the HM model for l =",l,"components","\n")
            estHMt <- lmestCont(responsesFormula = responsesFormula,
                                index = c(id, "time"),
                                data = Yt,
                                k = l,
                                modBasic = 1,
                                tol = 10^-10) 
            lkHMc[t] = estHMt$lk
            out <-long2matrices(id=Yt[,names(Yt)==id],time =Yt$time,Y =Yt[,names(Yt)%in%responses])
            cat("Estimate the constrained FM model for k =",l^3,"components","\n")
            estFMt <- estFM2(aperm(out$YY,c(1,3,2)),l=l, Mu0=t(estHMt$Mu),Si0=estHMt$Si,la0=estHMt$piv,Rho0=estHMt$Pi)
            lkFMc[t] = estFMt$lk
          }
          devc = 2*(lkFMc-lkHMc)
          pvc = 1-pchisq(devc,l^3-l^2)
          
          Pvc[,l] = pvc
          dft[l] = l^3-l^2
          df[l] = l^TT-l^2
          pvalueB[l] = min(1,min(pvc)*(TT-2))
          rejectB[l] = any(pvc<al)
        
          pvc1 <- sort(pvc)
          tmp <- (1:(TT-2))/(TT-2)
          pvalueS[l] <- min(1,min(pvc1/(1:(TT-2)))*(TT-2))
          rejectS[l] = any(pvc1<=tmp*al)
        }
      }
    }
  }

# ---- output ---
  lkHM = lkHM[lv]; names(lkHM) = lv
  aicHM = aicHM[lv]; names(aicHM) = lv
  bicHM = bicHM[lv]; names(bicHM) = lv
  if(type=="single"){
    Dev1=Dev1[lv]; names(Dev1) = lv
    pp = pvalue[lv]; names(pp)=lv
    reject = reject[lv]; names(reject) = lv  
    df = df[lv]; names(df) = lv  
    if(all(reject)){  
      selstates = max(lv)
    }else{  
      selstates = lv[which.min(reject==FALSE)]  
    }  
  }else{
    Dev1 = NULL
    ppB = pvalueB[lv]; names(ppB)=lv
    rejectB = rejectB[lv]; names(rejectB) = lv 
    ppS = pvalueS[lv]; names(ppS)=lv
    rejectS = rejectS[lv]; names(rejectS) = lv 
    df = df[lv]; dft = dft[lv]; names(df) = lv  
    if(all(rejectB)){  
      selstatesB = max(lv)
    }else{  
      selstatesB = lv[which.min(rejectB==FALSE)]  
    } 
    if(all(rejectS)){  
      selstatesS = max(lv)
    }else{  
      selstatesS = lv[which.min(rejectS==FALSE)]  
    } 
  }
  if(type=="single") out=list(lkHM=lkHM,aicHM=aicHM,bicHM=bicHM,Dev1=Dev1,n=n,TT=TT,pvalue=pp,reject= reject,
           df=df,selstates = selstates)   
  else if(TT>3 & type=="triplets") out=list(lkHM=lkHM,aicHM=aicHM,bicHM=bicHM,Dev1=Dev1,n=n,TT=TT,pvalueB=ppB,pvalueS=ppS,rejectB= rejectB,rejectS=rejectS,
                                            df=df,dft = dft,selstatesB = selstatesB,selstatesS=selstatesS,Pvc = Pvc[,lv,drop=FALSE])
  out$call = match.call()
  class(out) = "HM_FM_Test"
  return(out)

}