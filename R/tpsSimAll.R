tpsSimAll <-
function(B=1000,betaTruth,X,N,nII,interaction=NULL, ccDesign=NULL,
                      alpha=.05,threshold=c(-Inf,Inf),digits=NULL,
                      betaNames=NULL, referent=2, monitor=NULL,
                      cohort=TRUE,NI=NULL){
  beta <- betaTruth
  ##possible errors
  if(min(N,nII)<0){
    print("Error: sample size is not positive")
    return(-1)
  }

  ## case-control option
  if(cohort!=TRUE){
    if(length(NI)!=2){
      print("Error: 'NI' is a pair of Phase I sample sizes for controls and cases")
      return(-1)
    }
    if(min(NI)<0){
      print("Error: sample size is not positive")
      return(-1)
    }
    if(is.null(NI[1])){
      print("Warning: Phase I case-control sample size is not specified")
      cohort <- TRUE
    }
  } 

  if(min(X)<0)print("Error: invalid design matrix 'X'")
  if(length(nII)!=2){
    print("Error: 'nII' is not a pair of two numbers")
    return(-1)
  }
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

  
  if(length(interaction)>0){
    if(length(interaction)>=(ncol(X)-1)){
      print("Error: too many elements of 'interation'")
      return(-1)
    }
    if(length(setdiff(interaction,(2:ncol(X))))!=0){
      print("Error: invalid 'interaction'")
      return(-1)
    }
  }
  if(!is.null(threshold[1])){
    if(length(threshold)!=2){
      print("Error: 'threshold' is not a pair of numbers")
    }
  }


  ## check betaNames
  if(!is.null(threshold[1])){
    if(length(threshold)!=2){
      print("Error: 'threshold' is not a pair of numbers")
    }
  }
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

  ##check ccDesign
  if(!is.null(ccDesign)[1]){
    if(min(ccDesign)<0){
      print("Error: sample size is not positive")
      return(-1)
    }
    if(length(ccDesign)!=2){
      print("Error: invalid dimension of 'ccDesign'")
      return(-1)
    }
  }
  if(is.null(ccDesign[1])==1){
    ccDesign <- nII
  }

 
  # referent
  if(length(referent)!=1){
    print("Error: 'referent' is not a number'")
    return(-1)
  }else{
    if(!is.element(referent,1:2)){
        print("Error: invalid 'referent'")
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

  
  ## # of covariates
  num.cov <- ncol(X)-1
  beta.0.ind <- which(beta==0)

  if(num.cov==1){ #only one covariate
    num.grp <- length(unique(X[,2]))
    n0 <- round(rep(nII[1],times=num.grp)/num.grp)
    n1 <- round(rep(nII[2],times=num.grp)/num.grp)
    ## if there is one covariate, then strata = 2
    return(tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=2,n0=n0,n1=n1,ccDesign=ccDesign,alpha=alpha,threshold=threshold,digits=digits,referent=referent,betaNames=betaNames,monitor=monitor))
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
  
  
  dummy <- "mat.strata <- expand.grid(c(0,1)"
  if(num.cov > 1){
    for(i in 2:num.cov) dummy <- paste(dummy, " ,c(0,1)", sep="")
  }
  dummy <- paste(dummy, ")", sep="")
  cat(dummy, file="delete_this.R")
  source("delete_this.R")
  if(length(interaction)>0){
    mat.strata[,(interaction-1)] <- 0
    del.row <- NULL
    for(i in 1:(nrow(mat.strata)-1)){ 
      for(j in (i+1):(nrow(mat.strata))){
        if(prod(mat.strata[i,]==mat.strata[j,])!=0){
          if(j==nrow(mat.strata)){ #keep last  
            del.row <- c(del.row,i)
          }else{ 
            del.row <- c(del.row,j)
          }
        }
      }
    }
    mat.strata <- mat.strata[-del.row,]
  }
  num.des <- nrow(mat.strata)-2

  strata.name <- NULL
  result <- NULL
  result1 <- NULL
  result2 <- NULL

  #1st design
  strata <- sort(which(mat.strata[(1+1),]==1))+1
  num.grp <- 1
  strata.name.i <- NULL
  num.strata <- length(strata)
  for(j in 1:num.strata){
    num.grp <- num.grp*length(unique(X[,strata[j]]))
    strata.name.i <- paste(as.character(strata.name.i),strata[j], sep="")
  }
  strata.name <- c(strata.name,strata.name.i)
  n0 <- round(rep(nII[1],times=num.grp)/num.grp)
  n1 <- round(rep(nII[2],times=num.grp)/num.grp)
  print(paste("There are",(num.des+1),"two-phase designs to simulate."))
  result1 <- tpsSimAllOne(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,cohort=cohort,NI=NI)
  print("Design 1 complete.")
  
  #2nd design
  strata <- sort(which(mat.strata[(2+1),]==1))+1
  num.grp <- 1
  strata.name.i <- NULL
  num.strata <- length(strata)
  for(j in 1:num.strata){
    num.grp <- num.grp*length(unique(X[,strata[j]]))
    strata.name.i <- paste(as.character(strata.name.i),strata[j], sep="")
  }
  strata.name <- c(strata.name,strata.name.i)
  n0 <- round(rep(nII[1],times=num.grp)/num.grp)
  n1 <- round(rep(nII[2],times=num.grp)/num.grp)
  result2 <- tpsSimAllOne(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,cohort=cohort,NI=NI)
  print("Design 2 complete.")
  result <- list(result1,result2)
  
  if(num.cov > 2){
    for(i in 3:num.des){
      strata <- sort(which(mat.strata[(i+1),]==1))+1
      num.grp <- 1
      strata.name.i <- NULL
      num.strata <- length(strata)
      for(j in 1:num.strata){
        num.grp <- num.grp*length(unique(X[,strata[j]]))
        strata.name.i <- paste(as.character(strata.name.i),strata[j], sep="")
      }
      strata.name <- c(strata.name,strata.name.i)
      n0 <- round(rep(nII[1],times=num.grp)/num.grp)
      n1 <- round(rep(nII[2],times=num.grp)/num.grp)
      result[[i]] <- tpsSimAllOne(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,monitor=monitor,cohort=cohort,NI=NI)
      print(paste("Design",i,"complete."))
    }
  }
  result.cc <- NULL
  result.cc <- ccSimOne(B=B,betaTruth=beta,X=X,N=N,n0=ccDesign[1],n1=ccDesign[2],alpha=alpha,threshold=threshold,monitor=monitor,NI=NI)
  print(paste("Design",(i+1),"complete."))

  
  ##mean
  colname <- betaNames 
  rowname <- c("CD","CC",rep(c(" WL"," PL"," ML"),times=num.des))

  max.nchar <- max(nchar(as.numeric(strata.name)))
  spc <- NULL
  for(i in 1:max.nchar){
    spc <- paste(spc," ", sep="")
  }
  for(i in 1:num.des){
    for(j in 1:(max.nchar-1)){
      if(nchar(strata.name[i])==j){
        for(k in 1:(max.nchar-j)){
          strata.name[i] <- paste(strata.name[i]," ",sep="")
        }
        break
      }
    }
  }   
  j <- 1
  for(i in 1:(num.des*3)){
    if(i%%3==1){
      rowname[(i+2)] <- paste(strata.name[j],rowname[(i+2)],sep="")
      j <- j+1
    }else{
      rowname[(i+2)] <- paste(spc,rowname[(i+2)],sep="")
    }
  }
  
  matrix.mean <- NULL
  matrix.mean <- rbind(result[[1]]$mean[1,],result.cc$mean[2,],result[[1]]$mean[2:4,])
  for(i in 2:num.des){
    matrix.mean <- rbind(matrix.mean, result[[i]]$mean[-1,])
  }
  colnames(matrix.mean) <- colname
  rownames(matrix.mean) <- rowname


  ##bias in mean
  matrix.bias <- NULL
  matrix.bias <- rbind(result[[1]]$bias.mean[1,],result.cc$bias.mean[2,],result[[1]]$bias.mean[2:4,])
  for(i in 2:num.des){
    matrix.bias <- rbind(matrix.bias,result[[i]]$bias.mean[-1,])
  }
  rownames(matrix.bias) <- rowname
  colnames(matrix.bias) <- colname
  
  ##percent bias in mean
  matrix.pct.bias <- NULL
  matrix.pct.bias <- rbind(result[[1]]$pct.bias.mean[1,],result.cc$pct.bias.mean[2,],result[[1]]$pct.bias.mean[2:4,])
  for(i in 2:num.des){
    matrix.pct.bias <- rbind(matrix.pct.bias,result[[i]]$pct.bias.mean[-1,])
  }
  rownames(matrix.pct.bias) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.pct.bias) <- colname
  }else{
    colnames(matrix.pct.bias) <- colname[-beta.0.ind]
  }
  
  ##median
  matrix.med <- NULL
  matrix.med <- rbind(result[[1]]$median[1,],result.cc$median[2,],result[[1]]$median[2:4,])
  for(i in 2:num.des){
    matrix.med <- rbind(matrix.med,result[[i]]$median[-1,])
  }
  colnames(matrix.med) <- colname
  rownames(matrix.med) <- rowname

  ## bias in median
  matrix.bias.med <- NULL
  matrix.bias.med <- rbind(result[[1]]$bias.median[1,],result.cc$bias.median[2,],result[[1]]$bias.median[2:4,])
  for(i in 2:num.des){
    matrix.bias.med <- rbind(matrix.bias.med,result[[i]]$bias.median[-1,])
  }
  colnames(matrix.bias.med) <- colname
  rownames(matrix.bias.med) <- rowname
  
  
  ## Percent Bias in median
  matrix.pct.bias.med <- NULL
  matrix.pct.bias.med <- rbind(result[[1]]$pct.bias.med[1,],result.cc$pct.bias.med[2,],result[[1]]$pct.bias.med[2:4,])
  for(i in 2:num.des){
    matrix.pct.bias.med <- rbind(matrix.pct.bias.med,result[[i]]$pct.bias.med[-1,])
  }
  rownames(matrix.pct.bias.med) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.pct.bias.med) <- colname
  }else{
    colnames(matrix.pct.bias.med) <- colname[-beta.0.ind]
  }
  
  ## Standard Deviation
  matrix.sd <- NULL
  matrix.sd <- rbind(result[[1]]$sd[1,],result.cc$sd[2,],result[[1]]$sd[2:4,])
  for(i in 2:num.des){
    matrix.sd <- rbind(matrix.sd,result[[i]]$sd[-1,])
  }
  colnames(matrix.sd) <- colname
  rownames(matrix.sd) <- rowname
  
  ## Relative Uncertainty
  sd.comp <- matrix(rep(matrix.sd[referent,],times=nrow(matrix.sd)),byrow=TRUE,nc=length(beta),nr=nrow(matrix.sd))
  matrix.rel.unc <- matrix.sd/sd.comp*100
  matrix.rel.unc[referent,] <- 100.000
  colnames(matrix.rel.unc) <- colname
  rownames(matrix.rel.unc) <- rowname

  ## MSE
  matrix.mse <- NULL
  matrix.mse <- rbind(result[[1]]$mean.squared.error[1,],result.cc$mean.squared.error[2,],result[[1]]$mean.squared.error[2:4,])
  for(i in 2:num.des){
    matrix.mse <- rbind(matrix.mse,result[[i]]$mean.squared.error[-1,])
  } 
  colnames(matrix.mse) <- colname
  rownames(matrix.mse) <- rowname
  
  ## Mean Reperted Standard Error
  matrix.rep.sd <- NULL
  matrix.rep.sd <- rbind(result[[1]]$reported.standard.error[1,],result.cc$reported.standard.error[2,],result[[1]]$reported.standard.error[2:4,])
  for(i in 2:num.des){
    matrix.rep.sd <- rbind(matrix.rep.sd,result[[i]]$reported.standard.error[-1,])
  } 
  colnames(matrix.rep.sd) <- colname
  rownames(matrix.rep.sd) <- rowname
  
  ## Ratio of Reported and Actual Standard Deviation
  matrix.sd.rep.true <- NULL
  matrix.sd.rep.true <- rbind(result[[1]]$sd.reported.vs.actual[1,],result.cc$sd.reported.vs.actual[2,],result[[1]]$sd.reported.vs.actual[2:4,])
  for(i in 2:num.des){
    matrix.sd.rep.true <- rbind(matrix.sd.rep.true,result[[i]]$sd.reported.vs.actual[-1,])
  } 
  colnames(matrix.sd.rep.true) <- colname
  rownames(matrix.sd.rep.true) <- rowname
  
  ## Coverage Probability
  matrix.cov.prob <- NULL
  matrix.cov.prob <- rbind(result[[1]]$coverage.probability[1,],result.cc$coverage.probability[2,],result[[1]]$coverage.probability[2:4,])
  for(i in 2:num.des){
    matrix.cov.prob <- rbind(matrix.cov.prob,result[[i]]$coverage.probability[-1,])
  } 
  colnames(matrix.cov.prob) <- colname
  rownames(matrix.cov.prob) <- rowname

  ## Power
  matrix.power <- NULL
  matrix.power <- rbind(result[[1]]$power[1,],result.cc$power[2,],result[[1]]$power[2:4,])
  for(i in 2:num.des){
    matrix.power <- rbind(matrix.power,result[[i]]$power[-1,])
  } 
  rownames(matrix.power) <- rowname
  if(length(beta.0.ind)==0){
    colnames(matrix.power) <- colname
  }else{
    colnames(matrix.power) <- colname[-beta.0.ind]
  }
  
  ##NAs
  matrix.na <- NULL
  matrix.na <- cbind(result[[1]]$na[,1],result.cc$na[,2],result[[1]]$na[,2:4])
  for(i in 2:num.des){
    matrix.na <- cbind(matrix.na,result[[i]]$na[,-1])
  } 
  matrix.na <- t(matrix.na)
  colnames(matrix.na) <- c("Point Estimate","Above Threshold","Standard Deviation")
  rownames(matrix.na) <- rowname

  output <- NULL
  output$B <- B
  output$beta <- beta
  output$X <- X
  output$N <- N
  output$strata <- 0
  output$strata.des <- strata
  output$nII <- nII
  output$ccDesign <- ccDesign
  output$referent <- referent
  output$alpha <- alpha
  output$threshold <- threshold
  if(is.null(digits)){
    output$digits <- 0
  }else{
    output$digits <- digits
  }
  output$cohort <- cohort
  if(cohort==TRUE){
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

  class(output) <- "tpsSim"
    
  return(output)
}

