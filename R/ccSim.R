ccSim <-
function(B=1000,betaTruth,X,N,nCC,r,refDesign=1,
                  alpha=.05,threshold=c(-Inf,Inf),digits=NULL,
                  betaNames=NULL,monitor=NULL, NI=NULL){
  
  beta <- betaTruth
  if(length(refDesign)!=1){
    print("Error: 'refDesign is not a single number")
    return(-1)
  }
  if(refDesign<0){
    print("Error: 'refDesign' is not positive")
    return(-1)
  }

  if(length(r)==1){
    n0 <- round(nCC*r/(r+1))
    n1 <- nCC-n0
    return(tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=1,n0=n0,n1=n1,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor,NI=NI))
  }
  
  if(length(nCC)!=1){
    print("Error: 'nCC' is not a single number")
    return(-1)
  }
  if(min(N,nCC)<0){
    print("Error: sample size is not positive")
    return(-1)
  }

  if(min(X)<0)print("Error: invalid design matrix 'X'")
  if(min(r)<0)print("Error: invalid control-case ratio")

  
  if(length(N)!=nrow(X)){
    print("Error: invalid dimensions of 'N' and 'X'")
    return(-1)
  }

  ## valid beta
  nlevs <- 1
  for(i in 1 : (ncol(X)-1) ){
    nlevs <- c(nlevs,length(unique(X[ , (i+1)] ))-1)
  }
  if(length(beta)!=sum(nlevs)){
    print("Error: invalid dimension of 'beta'")
    return(-1)
  }  

  
  if(!is.null(threshold[1])){
    if(length(threshold)!=2){
      print("Error: 'threshold' is not a pair of numbers")
    }
  }
  

  ## check betaNames
  if(!is.null(betaNames[1])){
    if(length(beta)!=length(betaNames)){
      print("Error: 'beta' and 'betaNames' are not of the same length")
      return(-1)
    }
    if(!is.character(betaNames)){
      print("Error: elements of 'betaNames' are not character")
      return(-1)
    }
  }

  ## variable names, beta names
  if(!is.null(colnames(X))){
    var.name <- colnames(X)[-1]
    if(is.null(betaNames)){
      if(is.null(names(beta))){
        betaNames <- "Intercept"
        for(i in 2:ncol(X)){
          betaNames <- c(betaNames,paste("(",var.name[i-1],".",1:nlevs[i],")", sep=""))
        }
        names(beta) <- betaNames
      }else{
        betaNames <- names(beta)
      }
    }
  }else{
    colnames(X) <- c("Intercept",paste("V", 1:(ncol(X)-1), sep=""))
    var.name <- colnames(X)[-1]
    if(is.null(betaNames)){
      if(is.null(names(beta))){
        betaNames <- "Intercept"
        var.name2 <- paste("(V", 1:(ncol(X)-1),")", sep="")
        for(i in 2:ncol(X)){
          betaNames <- c(betaNames,paste(var.name[i-1],".",1:nlevs[i], sep=""))
        }
        names(beta) <- betaNames
      }else{
        betaNames <- names(beta)
      }
    }
  }  

  ##check monitor
  if(is.null(monitor)){
    monitor <- B+1
  }else{
    monitor <- as.integer(monitor)
    if(monitor<0){
      print("Error: 'monitor' is not a positive number")
      return(-1)
    }
  }

  ## NI
  if(!is.null(NI[1])){
    if(length(NI)!=2){
      print("Error: 'NI' is a pair of Phase I sample sizes for controls and cases")
      return(-1)
    }
    if(min(NI)<0){
      print("Error: sample size is not positive")
      return(-1)
    }
  }else{
    cohort <- TRUE
  }


  r <- c(refDesign,sort(unique(setdiff(r,refDesign))))
  num.des <- length(r)
  beta.0.ind <- which(beta==0)
  
  if(num.des==1){ #only one simulation
    num.grp <- length(unique(X[,2]))
    n0 <- round(nCC*r/(r+1))
    n1 <- nCC-n0
    return(tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=1,n0=n0,n1=n1,ccDesign=NULL,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor))
  }
  if(length(r)==1){
    return(tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=1,n0=n0,n1=n1,ccDesign=NULL,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor))
  }
  
  result <- NULL
  result1 <- NULL
  result2 <- NULL

  print(paste("There are",num.des,"case-control designs to simulate."))

  #1st design
  n0 <- round(nCC*r[1]/(r[1]+1))
  n1 <- nCC-n0
  result1 <- ccSimOne(B=B,betaTruth=beta,X=X,N=N,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,NI=NI)
  print("Design 1 complete.")

  #2nd design
  n0 <- round(nCC*r[2]/(r[2]+1))
  n1 <- nCC-n0
  result2 <- ccSimOne(B=B,betaTruth=beta,X=X,N=N,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,NI=NI)
  print("Design 2 complete.")
  result <- list(result1,result2)
  
  if(num.des > 2){
    for(i in 3:num.des){
      n0 <- round(nCC*r[i]/(r[i]+1))
      n1 <- nCC-n0
      result[[i]] <- ccSimOne(B=B,betaTruth=beta,X=X,N=N,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,NI=NI)
      print(paste("Design",i,"complete."))
    }
  }
 
  ##mean
  colname <- betaNames
  rowname <- c("CD",paste("CC:",r[1],"(ref)"),paste("CC:",r[-1]))
  matrix.mean <- NULL
  matrix.mean <- result[[1]]$mean
  for(i in 2:num.des){
    matrix.mean <- rbind(matrix.mean, result[[i]]$mean[2,])
  }
  colnames(matrix.mean) <- colname
  rownames(matrix.mean) <- rowname

  ##bias in mean
  matrix.bias <- NULL
  matrix.bias <- result[[1]]$bias.mean
  for(i in 2:num.des){
    matrix.bias <- rbind(matrix.bias,result[[i]]$bias.mean[2,])
  }
  rownames(matrix.bias) <- rowname
  colnames(matrix.bias) <- colname

  ##percent bias in mean
  matrix.pct.bias <- NULL
  matrix.pct.bias <- result[[1]]$pct.bias.mean
  for(i in 2:num.des){
    matrix.pct.bias <- rbind(matrix.pct.bias,result[[i]]$pct.bias.mean[2,])
  }
  rownames(matrix.pct.bias) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.pct.bias) <- colname
  }else{
    colnames(matrix.pct.bias) <- colname[-beta.0.ind]
  }
  
  ##median
  matrix.med <- NULL
  matrix.med <- result[[1]]$median
  for(i in 2:num.des){
    matrix.med <- rbind(matrix.med,result[[i]]$median[2,])
  }
  colnames(matrix.med) <- colname
  rownames(matrix.med) <- rowname

  ##bias in median
  matrix.bias.med <- NULL
  matrix.bias.med <- result[[1]]$bias.median
  for(i in 2:num.des){
    matrix.bias.med <- rbind(matrix.bias.med,result[[i]]$bias.median[2,])
  }
  colnames(matrix.bias.med) <- colname
  rownames(matrix.bias.med) <- rowname

  
  ## Percent Bias in median
  matrix.pct.bias.med <- NULL
  matrix.pct.bias.med <- result[[1]]$pct.bias.med
  for(i in 2:num.des){
    matrix.pct.bias.med <- rbind(matrix.pct.bias.med,result[[i]]$pct.bias.med[2,])
  }
  rownames(matrix.pct.bias.med) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.pct.bias.med) <- colname
  }else{
    colnames(matrix.pct.bias.med) <- colname[-beta.0.ind]
  }
  
  ## Standard Deviation
  matrix.sd <- NULL
  matrix.sd <- result[[1]]$sd
  for(i in 2:num.des){
    matrix.sd <- rbind(matrix.sd,result[[i]]$sd[2,])
  }
  colnames(matrix.sd) <- colname
  rownames(matrix.sd) <- rowname
  
  ## Relative Uncertainty
  matrix.rel.unc <- matrix.sd
  matrix.rel.unc <- matrix.rel.unc/matrix(rep(matrix.rel.unc[2,],times=nrow(matrix.sd)),byrow=TRUE,nr=nrow(matrix.sd))*100

  ## MSE
  matrix.mse <- NULL
  matrix.mse <- result[[1]]$mean.squared.error
  for(i in 2:num.des){
    matrix.mse <- rbind(matrix.mse,result[[i]]$mean.squared.error[2,])
  } 
  colnames(matrix.mse) <- colname
  rownames(matrix.mse) <- rowname
  
  ## Mean Reperted Standard Error
  matrix.rep.sd <- NULL
  matrix.rep.sd <- result[[1]]$reported.standard.error
  for(i in 2:num.des){
    matrix.rep.sd <- rbind(matrix.rep.sd,result[[i]]$reported.standard.error[2,])
  } 
  colnames(matrix.rep.sd) <- colname
  rownames(matrix.rep.sd) <- rowname
  
  ## Ratio of Reported and Actual Standard Deviation
  matrix.sd.rep.true <- NULL
  matrix.sd.rep.true <- result[[1]]$sd.reported.vs.actual
  for(i in 2:num.des){
    matrix.sd.rep.true <- rbind(matrix.sd.rep.true,result[[i]]$sd.reported.vs.actual[2,])
  } 
  colnames(matrix.sd.rep.true) <- colname
  rownames(matrix.sd.rep.true) <- rowname
  
  ## Coverage Probability
  matrix.cov.prob <- result[[1]]$coverage.probability
  for(i in 2:num.des){
    matrix.cov.prob <- rbind(matrix.cov.prob,result[[i]]$coverage.probability[2,])
  } 
  colnames(matrix.cov.prob) <- colname
  rownames(matrix.cov.prob) <- rowname

  ## Power
  matrix.power <- NULL
  matrix.power <- result[[1]]$power
  for(i in 2:num.des){
    matrix.power <- rbind(matrix.power,result[[i]]$power[2,])
  } 
  rownames(matrix.power) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.power) <- colname
  }else{
    colnames(matrix.power) <- colname[-beta.0.ind]
  }
  
  ##NAs
  matrix.na <- NULL
  matrix.na <- result[[1]]$na
  for(i in 2:num.des){
    matrix.na <- cbind(matrix.na,result[[i]]$na[,2])
  } 
  matrix.na <- t(matrix.na)
  colnames(matrix.na) <- c("Point Estimate","Above Threshold","Standard Deviation")
  rownames(matrix.na) <- rowname
  
  output <- NULL
  output$B <- B
  output$beta <- beta
  output$X <- X
  output$N <- N
  output$r <- r
  output$nCC <- nCC
  output$refDesign <- refDesign
  output$alpha <- alpha
  output$threshold <- threshold
  if(is.null(digits)){
    output$digits <- 0
  }else{
    output$digits <- digits
  }
  if(!is.null(NI[1])){
    output$NI <- NI
  }
  output$mean <- matrix.mean
  output$bias.mean <- matrix.bias
  output$pct.bias.mean <- matrix.pct.bias
  output$median <- matrix.med
  output$bias.median <- matrix.bias.med
  output$pct.bias.med <- matrix.pct.bias.med
  output$sd <- matrix.sd
  output$relative.uncertainty <- matrix.rel.unc 
  output$mean.squared.error <- matrix.mse
  output$reported.standard.error <- matrix.rep.sd
  output$sd.reported.vs.actual <- matrix.sd.rep.true
  output$coverage.probability <- matrix.cov.prob
  output$power <- matrix.power
  output$na <- matrix.na

  class(output) <- "ccSim"
    
  return(output)
}

