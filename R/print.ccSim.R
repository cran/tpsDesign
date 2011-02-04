print.ccSim <-
function(x,...,digits=-1){
  if(digits==-1){
    if(x$digits==0)digits <- 2
    else digits <- x$digits
  }
  cat("Number of simulations, B:",x$B,"\n")
  if(is.null(x$NI[1])){
    cat("Sample size at Phase I:",sum(x$N),"\n")
  }else{
    cat("Samle size at Phase I for controls:",x$NI[1],"\n")
    cat("Samle size at Phase I for cases:",x$NI[2],"\n")
  }
  cat("Sample size at Phase II, nCC:",x$nCC,"\n")
  cat("'True' regression coefficients, betaTruth:",x$beta,"\n")
  cat("\n")

  cat("Percent Bias in Mean\n")
  print(round(x$pct.bias.mean,digits=digits))
  cat("\n") 

  cat(1-x$alpha,"Level Coverage Probability\n")
  print(round(x$coverage.probability,digits=digits))
  cat("\n")

  cat("Relative Uncertainty\n")
  y <- x$relative.uncertainty
  y[,1] <- NA
  print(round(y,digits=digits))
  cat("\n") 
        
  cat("Number of NAs\n")
  if(x$threshold[2]!=Inf){
    print(x$na[-1,])
  }else{
    print(x$na[-1,-2])
  }
}

