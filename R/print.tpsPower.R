print.tpsPower <-
function(x,...,digits=-1){
  cat("Number of simulations, B:",x$B,"\n")
  cat("Sample size at Phase I:",sum(x$N),"\n")
  cat("Phase I stratification variable:",colnames(x$X)[x$strata],"\n")
  cat("'True' regression coefficients, betaTruth:",x$beta,"\n")
  cat("\n")
  
  if(digits==-1){
    if(x$digits==0)digits <- 2
    else digits <- x$digits
  }
    cat("Power\n")
  print(round(x$power,digits=digits))
  cat("\n") 

  cat("Number of NAs\n")
  if(x$threshold[2]!=Inf){
    print(x$na)
  }else{
    print(x$na[,-2])
  }
}
