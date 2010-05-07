tpsPower <-
function(B=1000,betaTruth,X,N,strata,nII,alpha=.05,threshold=c(-Inf,Inf),digits=NULL,betaNames=NULL,monitor=NULL){

  beta <- betaTruth
  ##possible errors
  if(min(N,nII)<0){
    print("Error: sample size is not positive")
    return(-1)
  }
  if(max(nII)>sum(N)){
    print("Error: one of 'nII' is more than sample size at phase I")
    return(-1)
  }
  if(min(X)<0)print("Error: invalid design matrix 'X'")
  
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
  

  if(length(intersect(strata,1:ncol(X)))==0){
    print("Error: 'strata' is invalid")
    return(-1)
  }
  
  if(length(strata)==1& strata[1]==1){
    print("Error: Use 'ccPower' for case control design")
    return(-1)
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

  result <- NULL
  result1 <- NULL
  result2 <- NULL
  
  nII <- sort(nII)
  num.nii <- length(nII)
  beta.0.ind <- which(beta==0)
  
  #1st design
  num.grp <- 1
  num.strata <- length(strata)
  for(j in 1:num.strata){
    num.grp <- num.grp*length(unique(X[,strata[j]]))
  }
  n0 <- round(rep(nII[1]/2,times=num.grp)/num.grp)
  n1 <- round(rep(nII[1]/2,times=num.grp)/num.grp)
  print(paste("There are",num.nii,"two-phase designs to simulate"))     
  result1 <- tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor)
  print("Design 1 complete")
  if(num.nii==1){
    output <- NULL
    output$B <- result1$B
    output$beta <- result1$beta
    output$X <- result1$X
    output$N <- result1$N
    output$strata <- result1$strata
    output$nII <- as.numeric(nII)
    output$alpha <- result1$alpha
    output$threshold <- result1$threshold
    if(is.null(digits)){
      output$digits <- 0
    }else{
      output$digits <- result1$digits
    }
    output$power <- result1$power
    output$na <- result1$na
    class(output) <- "tpsPower"
    return(output)
  }
  
  #2nd design
  n0 <- round(rep(nII[2]/2,times=num.grp)/num.grp)
  n1 <- round(rep(nII[2]/2,times=num.grp)/num.grp)
  result2 <- tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor)
  print("Design 2 complete")
  result <- list(result1,result2)
  
  if(num.nii > 2){
    for(i in 3:num.nii){
      n0 <- round(rep(nII[i]/2,times=num.grp)/num.grp)
      n1 <- round(rep(nII[i]/2,times=num.grp)/num.grp)
      result[[i]] <- tpsSim(B=B,betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,alpha=alpha,threshold=threshold,digits=digits,betaNames=betaNames,monitor=monitor)
      print(paste("Design",i,"complete."))
    }
  }

  colname <- betaNames
  rowname <- c("CD",rep(c(" CC"," WL"," PL"," ML"),times=num.nii))

  max.numnii <- max(nchar(nII))
  spc <- NULL
  for(i in 1:max.numnii){
    spc <- paste(spc," ", sep="")
  }
  for(i in 1:num.nii){
    for(j in 1:(max.numnii-1)){
      if(nchar(nII[i])==j){
        for(k in 1:(max.numnii-j)){
          nII[i] <- paste(nII[i]," ",sep="")
        }
        break
      }
    }
  }   
  j <- 1
  for(i in 1:(num.nii*4)){
    if(i%%4==1){
      rowname[(i+1)] <- paste(nII[j],rowname[(i+1)],sep="")
      j <- j+1
    }else{
      rowname[(i+1)] <- paste(spc,rowname[(i+1)],sep="")
    }
  }
  
  ## Power
  matrix.power <- NULL
  matrix.power <- result[[1]]$power
  for(i in 2:num.nii){
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
  matrix.na <- result[[1]]$na
  for(i in 2:num.nii){
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
  output$strata <- strata
  output$nII <- as.numeric(nII)
  output$alpha <- alpha
  output$threshold <- threshold
  if(is.null(digits)){
    output$digits <- 0
  }else{
    output$digits <- digits
  }
  output$power <- matrix.power
  output$na <- matrix.na
  
  class(output) <- "tpsPower"
  
  return(output)
}

