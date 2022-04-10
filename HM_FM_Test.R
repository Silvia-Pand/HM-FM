HM_FM_Test <- function(data=data,id = "id",time="time",lv=2:4,responses=NULL,type=c("single","triplets")){

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
#             "triplets" = the multiple test based on consecutive triplets of observations)
#
# OUTPUT:
#lkHM    = value of the  log-likelihood of the HM model for each l in lk
#aicHM   = value of AIC index for the HM model for each l in lv
#bicHM   = value of BIC index for the HM model for each l in lv
#n       = sample size
#TT      = number of time occasions
#pvalue  = pvalue of the standard LR test for each l in lv (if type="single")   
#or overall pvalue of the multiple LR test for each l in lv (if type="triplets")

# ---- preliminaries ----
  type <- match.arg(type, choices = eval(formals(HM_FM_Test)$type))
  if(is.null(responses)) responses = colnames(data[,3:dim(data)[2]])
  responsesFormula <- lmestFormula(data=data,response=responses)$responsesFormula
 
  aicHM <-bicHM <- lkHM <- pvalue <- rep(NA,max(lv))
  out <-long2matrices(id=data[,names(data)==id],time =data[,names(data)==time],Y =data[,names(data)%in%responses])
  n <- dim(out$YY)[1]
  TT <- dim(out$YY)[2]
  YY <- out$YY
  
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
      if(type=="single"){
        print("Perform the standard LR test")
        #fit constrained FM model  
        cat("Estimate the FM2 model for k =",l^TT,"components","\n")
        estFM <- estFM2(aperm(YY,c(1,3,2)),l=l, Mu=t(estHM$Mu),Si=estHM$Si,piv0=estHM$piv,Pi0=estHM$Pi)
        lkFM <- estFM$lk
        Dev1 <- -2*(lkHM[l] - lkFM)
        pvalue[l] <- 1-pchisq(Dev1,l^TT-l^2)   
        
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
            estFMt <- estFM2(aperm(out$YY,c(1,3,2)),l=l, Mu=t(estHMt$Mu),Si=estHMt$Si,piv0=estHMt$piv,Pi0=estHMt$Pi)
            lkFMc[t] = estFMt$lk
          }
          devc = 2*(lkFMc-lkHMc)
          pvc = 1-pchisq(devc,l^3-l^2)
          pvalue[l] = min(1,min(pvc)*(TT-2))
        }
          
      }
    }
  }

# ---- output ---
  pp = pvalue[lv]
  names(pp)=lv
  out=list(lkHM=lkHM,aicHM=aicHM,bicHM=bicHM,n=n,TT=TT,pvalue=pp)
  return(out)

}