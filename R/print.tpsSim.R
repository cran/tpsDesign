print.tpsSim <-
function(x,...,digits=-1,referent=NULL){
  if(length(x$strata)!=1 | (!is.element(0,x$strata)&!is.element(1,x$strata))){
    if(digits==-1){
      if(x$digits==0)digits <- 2
      else digits <- x$digits
    }
    if(x$cohort==TRUE){
      cat("Cohort sampling\n")
    }else{
      cat("Case-control sampling\n")
    }
    cat("Number of simulations:",x$B,"\n")
    cat("Phase I stratification variable:",colnames(x$X)[x$strata],"\n")
    if(x$cohort==TRUE){
      cat("Sample size at Phase I:",sum(x$N),"\n")
    }else{
      cat("Samle size at Phase I for controls:",x$NI[1],"\n")
      cat("Samle size at Phase I for cases:",x$NI[2],"\n")
    }
    cat("Sample size at Phase II for controls, n0:",x$n0,"\n")
    cat("Sample size at Phase II for cases, n1:",x$n1,"\n")
    cat("Sample sizes for a referent case-control design, ccDesign:",x$ccDesign,"\n")
    cat("'True' regression coefficients, betaTruth:",x$beta,"\n")
    cat("\n")
    
    cat("Percent Bias in Mean\n")
    print(round(x$pct.bias.mean,digits=digits))
    cat("\n") 

    cat(1-x$alpha,"Level Coverage Probability\n")
    print(round(x$coverage.probability,digits=digits))
    cat("\n")

    y <- x$relative.uncertainty
    if(is.null(referent)||referent==x$referent){
      cat("Relative Uncertainty\n")
      if(x$referent!=2){
        y[2,1] <- NA
        print(round(y,digits=digits))
      }else{
        y[,1] <- NA
        print(round(y,digits=digits))
      }
    }else{
      if(!is.element(referent,1:5)){
        print("Warning: invalid 'referent'")
        cat("Relative Uncertainty\n")
        if(x$referent!=2){
          y[2,1] <- NA
          print(round(y,digits=digits))
        }else{
          y[,1] <- NA
          print(round(y,digits=digits))
        }
      }else{
        dim.beta <- length(x$beta)
        sd.comp <- rbind(x$sd[referent,],x$sd[referent,],x$sd[referent,],x$sd[referent,],x$sd[referent,])
        
        rel.unc <- x$sd/sd.comp*100
        rel.unc[referent,] <- 100.000
        if(referent==2){
          rel.unc[,1] <- NA
        }else{
          rel.unc[2,1] <- NA
        }
        cat("Relative Uncertainty\n")
        print(round(rel.unc,digits=digits))

      }
    }
    
    cat("\n")    
    cat("Number of NAs\n")
    if(x$threshold[1]!=-Inf&x$threshold[2]!=Inf){
      print(x$na)
    }else{
      print(x$na[-2,])
    }
    
  }else if(length(x$strata)==1 && x$strata==1){ ##CC
    if(digits==-1){
      if(x$digits==0)digits <- 2
      else digits <- x$digits
    }
    
    cat("Number of simulations, B:",x$B,"\n")
    cat("Case-control:\n")
    if(is.null(x$NI[1])){
      cat("Sample size at Phase I:",sum(x$N),"\n")
    }else{
      cat("Samle size at Phase I for controls:",x$NI[1],"\n")
      cat("Samle size at Phase I for cases:",x$NI[2],"\n")
    }
    cat("Sample size for controls, n0:",x$n0,"\n")
    cat("Sample size for cases, n1:",x$n1,"\n")
    if(nrow(x$mean)==3){
      cat("Sample sizes for a referent case-control design, ccDesign:",x$ccDesign,"\n")
    }
    cat("'True' regression coefficients, betaTruth:",x$beta,"\n")
   
    cat("\n")
    
    cat("Percent Bias in Mean\n")
    print(round(x$pct.bias.mean,digits=digits))
    cat("\n") 

    cat(1-x$alpha,"Level Coverage Probability\n")
    print(round(x$coverage.probability,digits=digits))
    cat("\n") 

    y <- x$relative.uncertainty
    if(is.null(referent)||referent==x$referent){
      cat("Relative Uncertainty\n")
      if(x$referent==1){
        y[-1,1] <- NA
        print(round(y,digits=digits))
      }else{
        y[,1] <- NA
        print(round(y,digits=digits))
      }
    }else{
      if(!is.element(referent,1:3)){
        print("Warning: invalid 'referent'")
        cat("Relative Uncertainty\n")
        if(x$referent==1){
          y[-1,1] <- NA
          print(round(y,digits=digits))
        }else{
          y[,1] <- NA
          print(round(y,digits=digits))
        }
      }else{
        dim.beta <- length(x$beta)
        sd.comp <- rbind(x$sd[referent,],x$sd[referent,],x$sd[referent,])
        rel.unc <- x$sd/sd.comp*100
        rel.unc[referent,] <- 100.000
        if(referent==1){
          rel.unc[-1,1] <- NA
        }else{
          rel.unc[,1] <- NA
        }
        cat("Relative Uncertainty\n")        
        print(round(rel.unc,digits=digits))
      }
    }
    
    cat("\n")
    cat("Number of NAs\n")
    if(x$threshold[1]!=-Inf & x$threshold[2]!=Inf){
      print(x$na)
    }else{
      print(x$na[-2,])
    }
    
  }else if(length(x$strata)==1 && x$strata==0){ #all designs
    if(digits==-1){
      if(x$digits==0)digits <- 2
      else digits <- x$digits
    }
    cat("Number of simulations, B:",x$B,"\n")
    cat("Sample size at Phase I:",sum(x$N),"\n")
    cat("Sample size at Phase II for controls, n0:",x$nII[1],"\n")
    cat("Sample size at Phase II for cases, n1:",x$nII[2],"\n")
    cat("Sample size for a referent case-control design, ccDeisgn:",x$ccDesign,"\n")
    cat("'True' regression coefficients, betaTruth:",x$beta,"\n")
    cat("\n")

    colname <- paste(2:ncol(x$X),":=",colnames(x$X)[-1],"  ")
    cat("Percent Bias in Mean\n")
    cat(colname)
    cat("\n") 
    print(round(x$pct.bias.mean,digits=digits))
    cat("\n") 

    cat(1-x$alpha,"Level Coverage Probability\n")
    cat(colname)
    cat("\n") 
    print(round(x$coverage.probability,digits=digits))
    cat("\n")
    
    y <- x$relative.uncertainty
    if(is.null(referent)||referent==x$referent){
      cat("Relative Uncertainty\n")
      cat(colname)
      cat("\n")
      if(x$referent!=2){
        y[2,1] <- NA
        print(round(y,digits=digits))
      }else{
        y[,1] <- NA
        print(round(y,digits=digits))
      }
    }else{
      if(!is.element(referent,c(1,2))){
        print("Warning: invalid 'referent'")
        cat("Relative Uncertainty\n")
        cat(colname)
        cat("\n")
        if(x$referent!=2){
          y[2,1] <- NA
          print(round(y,digits=digits))
        }else{
          y[,1] <- NA
          print(round(y,digits=digits))
        }
      }else{
        dim.beta <- length(x$beta)
        nrows <- dim(x$sd)[1]
        sd.comp <- matrix(rep(x$sd[referent,],times=nrows),byrow=TRUE,nc=dim.beta,nr=nrows)
        rel.unc <- x$sd/sd.comp*100
        rel.unc[referent,] <- 100.000
        if(referent==2){
          rel.unc <- rel.unc[,1] <- NA
        }else{
          rel.unc[2,1] <- NA
        }
        cat("Relative Uncertainty\n")
        cat(colname)
        cat("\n") 
        print(round(rel.unc,digits=digits))
      }
    }

    cat("\n")
    cat("Number of NAs\n")
    cat(colname)
    cat("\n")
    if(x$threshold[1]!=-Inf & x$threshold[2]!=Inf){
      print(x$na)
    }else{
      print(x$na[,-2])
    }
  }
}

