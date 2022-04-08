sqg <- function(lev){

#        S = sqg(lev)
#
# generalization of sq; the maximum of the elements in
# any column of S is given by the corresponding entry of lev
# initial level is 0

  J = length(lev)
  lev = pmax(lev,0)
  if(J==1){
    S = matrix(0:lev[1],lev[1]+1,1)
  }else{
    Te = sqg(lev[2:J])
    nt = nrow(Te)
    S = cbind((0:lev[1])%x%rep(1,nt),rep(1,lev[1]+1)%x%Te)
  }
  S
  
}