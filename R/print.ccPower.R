print.ccPower <-
function(x,...,digits=-1){
  cat("Number of simulations, B: ",x$B,"\n")
  if(is.null(x$NI[1])){
    cat("Sample size at Phase I:",sum(x$N),"\n")
  }else{
    cat("Sample size at Phase I for controls:",x$NI[1],"\n")
    cat("Sample size at Phase I for casess:",x$NI[2],"\n")
  }
  cat("'True' regession coefficients, betaTruth:",x$beta,"\n")
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

