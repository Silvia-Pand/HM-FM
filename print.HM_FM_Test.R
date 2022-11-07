print.HM_FM_Test <-function(x, ...){
   cat("Call:\n")
   print(x$call)
   cat("\nTest restuls:\n")
   if(!is.null(x$pvalue)) print(cbind(df=x$df, "p-value"=x$pvalue))
   else print(cbind(df=x$df,"df (triplets)"= x$dft,"Bonferroni p-value"=x$pvalueB,"Simes p-value" = x$pvalueS))
 } 