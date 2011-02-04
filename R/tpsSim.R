tpsSim <-
function(B=1000, betaTruth, X, N, strata, n0=NULL, n1=NULL, nII=NULL, interaction=NULL,
                   ccDesign=NULL, alpha=.05, threshold=c(-Inf,Inf),
                   digits=NULL, betaNames=NULL,referent=2, monitor=NULL,
                   cohort=TRUE, NI=NULL){

  beta <- betaTruth
  ## valid strata
  if(length(intersect(strata,0:ncol(X)))==0){
    print("Error: 'strata' is invalid")
    return(-1)
  }
  ## all design or specific?
  if(is.element(0,strata)){
    if(length(strata)>1){
      print("Error: 'strata' is invalid")
      return(-1)
    }
    if(!is.null(n0)|!is.null(n1)){
      print("Warning: 'n0' and 'n1' are not used")
    }
    return(tpsSimAll(B=B,betaTruth=beta,X=X,N=N,nII=nII,interaction=interaction,
                      ccDesign=ccDesign,alpha=alpha,threshold=threshold,digits=digits,
                      referent=referent,monitor=monitor,cohort=cohort,NI=NI))
  }
  if(!is.null(nII[1]))print("Warning: 'nII' is not used")
  ## case-control or specific
  if(is.element(1,strata) & length(strata)>1){
    print("Error: 'strata' is invalid")
    return(-1)
  }


  ## valid N, n0 ,n1
  if(min(N,n0,n1)<0){
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
  if(!is.null(NI[1])&cohort==TRUE){
    print("Warning: 'cohort' option was chosen")
  }
         
  
  ##valid X
  if(min(X)<0)print("Error: invalid design matrix 'X'")

  ## valid ccDesign
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

  ## valid N and X   
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

  n.st <- 1
  for(i in 1:length(strata)){
    n.st <- n.st*length(unique(X[,strata[i]]))
  }
  if(length(n0)!=n.st){
    print("Error: invalid dimension of 'n0'")
    return(-1)
  }
  if(length(n1)!=n.st){
    print("Error: invalid dimension of 'n1'")
    return(-1)
  }
  rm(n.st)


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

  ##check referent
  if(length(referent)!=1){
    print("Error: 'referent' is not a number'")
    return(-1)
  }else{
    if(length(strata)==1&strata[1]==1){
      if(!is.element(referent,1:3)){
        print("Error: invalid 'referent'")
        return(-1)
      }
    }else if(!is.element(referent,1:5)){
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

  dim.beta <- length(beta)
  strata <- sort(strata)
    
  ## simulation
  if(strata[1]!=1){ #non CC
    ##simulation 
    result <- mat.or.vec(nr=dim.beta*8,nc=B)
    for(i in 1:B){
      if(i %% monitor==0)
        cat("Repetition",i,"of",B,"for a two-phase design complete\n") 
      result[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1,var.name=var.name,cohort=cohort,NI=NI)
    }
    
    ##simulation for case-control
    result.cc <- mat.or.vec(nr=dim.beta*8,nc=B)
    if(is.null(ccDesign[1])==1){ #not specified by user
      ccDesign <- c(sum(n0),sum(n1))
    }
    for(i in 1:B){
      if(i %% monitor==0)
        cat("Repetition",i,"of",B,"for a case-control design complete\n") 
      result.cc[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=1,n0=ccDesign[1],n1=ccDesign[2],cohort=cohort,NI=NI)
    }
    ##NA
    na <- numeric(12)
    na.cc <- numeric(3)
    
    ##NA by too large deviation of estimates from truth possibly due to numerical difficulty
    na[5] <- sum(apply(result[1:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
    na[6] <- sum(apply(result[(dim.beta+1):(2*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
    na[7] <- sum(apply(result[(2*dim.beta+1):(3*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
    na[8] <- sum(apply(result[(3*dim.beta+1):(4*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
    na.cc[2] <- sum(apply(result.cc[2:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta[-1],threshold=threshold)) #intercept estimate is not reliable in case-control
    
 
    ##NA by sigularity of design matrix
    na[9] <- sum(is.na(result[(4*dim.beta+1),]))
    na[10] <- sum(is.na(result[(5*dim.beta+1),]))
    na[11] <- sum(is.na(result[(6*dim.beta+1),]))
    for(i in 1:dim.beta){
      na[12] <- max(na[12],sum(is.na(result[(7*dim.beta+i),])))
    }
    na.cc[3] <- sum(is.na(result.cc[(dim.beta+1),]))

    ##NA for all estimates in each method
    na[1] <- sum(is.na(result[1,]))
    na[2] <- sum(is.na(result[(dim.beta+1),]))
    na[3] <- sum(is.na(result[(2*dim.beta+1),]))
    na[4] <- sum(is.na(result[(3*dim.beta+1),]))
    na.cc[1] <- sum(is.na(result.cc[2,]))
    
    ##NA removed
    result.original <- result
    result.cc.original <- result.cc
    
    for(i in 1:dim.beta){
      result[c(1:dim.beta,(4*dim.beta+1):(5*dim.beta)),
             result[i,]-beta[i] < threshold[1] |
             result[i,]-beta[i] > threshold[2] |
             is.na(result[(4*dim.beta+i),])
             | is.na(result[i,])] <- NA
      
      result[c((dim.beta+1):(2*dim.beta),(5*dim.beta+1):(6*dim.beta)),
             result[(dim.beta+i),]-beta[i] < threshold[1]
             |result[(dim.beta+i),]-beta[i] > threshold[2]
             | is.na(result[(5*dim.beta+i),])
             | is.na(result[(dim.beta+i),])] <- NA
      
      result[c((2*dim.beta+1):(3*dim.beta),(6*dim.beta+1):(7*dim.beta)),
             result[(2*dim.beta+i),]-beta[i] < threshold[1]
             |result[(2*dim.beta+i),]-beta[i] > threshold[2]
             | is.na(result[(6*dim.beta+i),]) | is.na(result[(2*dim.beta+i),])] <- NA
      
      result[c((3*dim.beta+1):(4*dim.beta),(7*dim.beta+1):(8*dim.beta)),
             result[(3*dim.beta+i),]-beta[i] < threshold[1]
             |result[(3*dim.beta+i),]-beta[i] > threshold[2]
             | is.na(result[(7*dim.beta+i),])
             | is.na(result[(3*dim.beta+i),])] <- NA

      if(i!=1){ #intercept estimate in case-control is not reliable
        result.cc[c(1:(2*dim.beta)),
                  result.cc[i,]-beta[i] < threshold[1]
                  |result.cc[i,]-beta[i] > threshold[2]
                  | is.na(result.cc[(dim.beta+i),])
                  | is.na(result.cc[i,])] <- NA
      }
    }

    ##mean
    est.ave <- apply(result[1:(4*dim.beta),],1,mean,na.rm=TRUE)
    est.ave.cc <- apply(result.cc[1:dim.beta,],1,mean,na.rm=TRUE)
  
    ##percent bias in mean
    matrix.beta <- rbind(beta,beta,beta,beta)
    matrix.bias.mean <- matrix(est.ave,nrow=4,byrow=TRUE)-matrix.beta
    matrix.bias.mean.cc <- est.ave.cc-beta
    pct.bias <- matrix.bias.mean/abs(matrix.beta)*100
    pct.bias.cc <- matrix.bias.mean.cc/abs(beta)*100
    beta.0.ind <- which(beta==0) #consideration of case of dividing by zero
    if(!length(beta.0.ind)){
      num.null.0 <- dim.beta
    }else{  
      num.null.0 <- dim.beta - length(beta.0.ind)
    }
    if(length(beta.0.ind) > 0){
      pct.bias <- pct.bias[,-beta.0.ind]
      pct.bias.cc <- pct.bias.cc[-beta.0.ind]
    }
 
    ##median
    med <- apply(result[1:(4*dim.beta),],1,median,na.rm=TRUE)
    med.cc <- apply(result.cc[1:dim.beta,],1,median,na.rm=TRUE)

    ##percent bias in median
    matrix.diff <- matrix(med,nrow=4,byrow=TRUE)-matrix.beta
    matrix.diff.cc <- med.cc-beta
    pct.bias.med <- matrix.diff/abs(matrix.beta)*100
    pct.bias.med.cc <- matrix.diff.cc/abs(beta)*100
    if(length(beta.0.ind) > 0){
      pct.bias.med <- pct.bias.med[,-beta.0.ind]
      pct.bias.med.cc <- pct.bias.med.cc[-beta.0.ind]
    }    
    
    ##standard deviation
    sds <- apply(result[1:(4*dim.beta),],1,sd,na.rm=TRUE)
    sds.cc <- apply(result.cc[1:dim.beta,],1,sd,na.rm=TRUE)
    
    ##relative uncertainty
    rel.unc <- rbind(sds[1:dim.beta],sds.cc,sds[(dim.beta+1):(2*dim.beta)],sds[(2*dim.beta+1):(3*dim.beta)],sds[(3*dim.beta+1):(4*dim.beta)])
    rel.unc <- rel.unc/matrix(rep(rel.unc[referent,],times=5),byrow=TRUE,nrow=5)*100
    rel.unc[referent,] <- 100.000
    
    ##sd average
    sd.ave <- apply(result[(4*dim.beta+1):(8*dim.beta),],1,mean,na.rm=TRUE)
    sd.ave.cc <- apply(result.cc[(dim.beta+1):(2*dim.beta),],1,mean,na.rm=TRUE)

    ##ratio of reported sds and the actual truth
    rel.sd <- sd.ave/sds
    rel.sd.cc <- sd.ave.cc/sds.cc
    
    ##MSE
    mse <- (est.ave-beta)^2+sds^2
    mse.cc <- (est.ave.cc-beta)^2+sds.cc^2
    
    ##coverage probability
    cov.prob <- 100*apply(ifelse(result[1:(4*dim.beta),] + qnorm(alpha/2)*result[(4*dim.beta+1):(8*dim.beta),] < beta &
                           beta < result[1:(4*dim.beta),]+qnorm(1-alpha/2)*result[(4*dim.beta+1):(8*dim.beta),]
                           ,1,0),1,mean,na.rm=TRUE)
  
    cov.prob.cc <- 100*apply(ifelse(result.cc[1:dim.beta,] + qnorm(alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),] < beta &
                           beta < result.cc[1:dim.beta,]+qnorm(1-alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                           ,1,0),1,mean,na.rm=TRUE)

    ##power with H0:beta=0
    power.0 <- 100-100*apply(ifelse(result[1:(4*dim.beta),] + qnorm(alpha/2)*result[(4*dim.beta+1):(8*dim.beta),] < 0 &
                              0 < result[1:(4*dim.beta),] + qnorm(1-alpha/2)*result[(4*dim.beta+1):(8*dim.beta),]
                              ,1,0),1,mean,na.rm=TRUE)
    power.0.cc <- 100-100*apply(ifelse(result.cc[1:dim.beta,]+ qnorm(alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),] < 0 &
                                 0 < result.cc[1:dim.beta,]+qnorm(1-alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                                 ,1,0),1,mean,na.rm=TRUE)
    if(length(beta.0.ind) > 0){
      power.0 <- power.0[-c(beta.0.ind,(beta.0.ind+dim.beta),(beta.0.ind+2*dim.beta),(beta.0.ind+3*dim.beta))]
      power.0.cc <- power.0.cc[-beta.0.ind]
    }
    

    cov.name <- betaNames
    
    matrix.mean <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.mean) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.mean) <- cov.name
    matrix.mean[1,] <- est.ave[1:dim.beta]
    matrix.mean[2,] <- est.ave.cc
    matrix.mean[3,] <- est.ave[(dim.beta+1):(2*dim.beta)]
    matrix.mean[4,] <- est.ave[(2*dim.beta+1):(3*dim.beta)]
    matrix.mean[5,] <- est.ave[(3*dim.beta+1):(4*dim.beta)]    
    matrix.mean[2,1] <- NA
    
    matrix.bias <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.bias) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.bias) <- cov.name
    matrix.bias[1,] <- matrix.bias.mean[1,]
    matrix.bias[2,] <- matrix.bias.mean.cc
    matrix.bias[3,] <- matrix.bias.mean[2,]
    matrix.bias[4,] <- matrix.bias.mean[3,]
    matrix.bias[5,] <- matrix.bias.mean[4,]
    matrix.bias[2,1] <- NA
    
    if(num.null.0!=0){
      matrix.pctbias <- as.matrix(mat.or.vec(nc=num.null.0,nr=5))
      rownames(matrix.pctbias) <- c("CD","CC","WL","PL","ML")
      if(!length(beta.0.ind)){
        colnames(matrix.pctbias) <- cov.name
      }else{
        colnames(matrix.pctbias) <- cov.name[-beta.0.ind]
      }
      matrix.pctbias[1,] <- pct.bias[1,]
      matrix.pctbias[2,] <- pct.bias.cc[1:num.null.0]  
      matrix.pctbias[3,] <- pct.bias[2,]
      matrix.pctbias[4,] <- pct.bias[3,]
      matrix.pctbias[5,] <- pct.bias[4,]
      matrix.pctbias[2,1] <- NA
    }

    matrix.med <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.med) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.med) <- cov.name
    matrix.med[1,] <- med[1:dim.beta]
    matrix.med[2,] <- med.cc
    matrix.med[3,] <- med[(dim.beta+1):(2*dim.beta)]
    matrix.med[4,] <- med[(2*dim.beta+1):(3*dim.beta)]
    matrix.med[5,] <- med[(3*dim.beta+1):(4*dim.beta)]    
    matrix.med[2,1] <- NA

    matrix.bias.med <- matrix.med
    matrix.beta <- rbind(beta,beta,beta,beta,beta)
    matrix.bias.med <- matrix.bias.med-matrix.beta
    
    if(num.null.0!=0){
      matrix.pctbias.med <- as.matrix(mat.or.vec(nc=num.null.0,nr=5))
      rownames(matrix.pctbias.med) <- c("CD","CC","WL","PL","ML")
      if(!length(beta.0.ind)){
        colnames(matrix.pctbias.med) <- cov.name
      }else{
        colnames(matrix.pctbias.med) <- cov.name[-beta.0.ind]
      }
      matrix.pctbias.med[1,] <- pct.bias.med[1,]
      matrix.pctbias.med[2,] <- pct.bias.med.cc[1:num.null.0]  
      matrix.pctbias.med[3,] <- pct.bias.med[2,]
      matrix.pctbias.med[4,] <- pct.bias.med[3,]
      matrix.pctbias.med[5,] <- pct.bias.med[4,]
      matrix.pctbias[2,1] <- NA
    }
    
    matrix.sd <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.sd) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.sd) <- cov.name
    matrix.sd[1,] <- sds[1:dim.beta]
    matrix.sd[2,] <- sds.cc
    matrix.sd[3,] <- sds[(dim.beta+1):(2*dim.beta)]
    matrix.sd[4,] <- sds[(2*dim.beta+1):(3*dim.beta)]
    matrix.sd[5,] <- sds[(3*dim.beta+1):(4*dim.beta)]
    matrix.sd[2,1] <- NA
    
    matrix.relunc <- rel.unc
    rownames(matrix.relunc) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.relunc) <- cov.name
    
    matrix.mse <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.mse) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.mse) <- cov.name
    matrix.mse[1,] <- mse[1:dim.beta]
    matrix.mse[2,] <- mse.cc
    matrix.mse[3,] <- mse[(dim.beta+1):(2*dim.beta)]
    matrix.mse[4,] <- mse[(2*dim.beta+1):(3*dim.beta)]
    matrix.mse[5,] <- mse[(3*dim.beta+1):(4*dim.beta)]
    matrix.mse[2,1] <- NA
    
    matrix.sdave <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.sdave) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.sdave) <- cov.name
    matrix.sdave[1,] <- sd.ave[1:dim.beta]
    matrix.sdave[2,] <- sd.ave.cc
    matrix.sdave[3,] <- sd.ave[(dim.beta+1):(2*dim.beta)]
    matrix.sdave[4,] <- sd.ave[(2*dim.beta+1):(3*dim.beta)]
    matrix.sdave[5,] <- sd.ave[(3*dim.beta+1):(4*dim.beta)]
    matrix.sdave[2,1] <- NA
    
    matrix.rel.sd <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.rel.sd) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.rel.sd) <- cov.name
    matrix.rel.sd[1,] <- rel.sd[1:dim.beta]
    matrix.rel.sd[2,] <- rel.sd.cc
    matrix.rel.sd[3,] <- rel.sd[(dim.beta+1):(2*dim.beta)]
    matrix.rel.sd[4,] <- rel.sd[(2*dim.beta+1):(3*dim.beta)]
    matrix.rel.sd[5,] <- rel.sd[(3*dim.beta+1):(4*dim.beta)]
    matrix.rel.sd[2,1] <- NA
    
    matrix.covprob <- mat.or.vec(nc=dim.beta,nr=5)
    rownames(matrix.covprob) <- c("CD","CC","WL","PL","ML")
    colnames(matrix.covprob) <- cov.name
    matrix.covprob[1,] <- cov.prob[1:dim.beta]
    matrix.covprob[2,] <- cov.prob.cc
    matrix.covprob[3,] <- cov.prob[(dim.beta+1):(2*dim.beta)]
    matrix.covprob[4,] <- cov.prob[(2*dim.beta+1):(3*dim.beta)]
    matrix.covprob[5,] <- cov.prob[(3*dim.beta+1):(4*dim.beta)]
    matrix.covprob[2,1] <- NA
    
    if(num.null.0!=0){
      matrix.power <- as.matrix(mat.or.vec(nc=num.null.0,nr=5))
      rownames(matrix.power) <- c("CD","CC","WL","PL","ML")
      if(!length(beta.0.ind)){
        colnames(matrix.power) <- cov.name
      }else{
        colnames(matrix.power) <- cov.name[-beta.0.ind]
      }
      matrix.power[1,] <- power.0[1:num.null.0]
      matrix.power[2,] <- power.0.cc[1:num.null.0]
      matrix.power[3,] <- power.0[(num.null.0+1):(2*num.null.0)]
      matrix.power[4,] <- power.0[(2*num.null.0+1):(3*num.null.0)]
      matrix.power[5,] <- power.0[(3*num.null.0+1):(4*num.null.0)]
      matrix.power[2,1] <- NA
    }
    
    matrix.na <- mat.or.vec(nc=5,nr=3)
    colnames(matrix.na) <- c("CD","CC","WL","PL","ML")
    rownames(matrix.na) <- c("Point Estimate","Above Threshold","Standard Deviation")
    matrix.na[1,] <- c(na[1],na.cc[1],na[2:4])
    matrix.na[2,] <- c(na[5],na.cc[2],na[6:8])
    matrix.na[3,] <- c(na[9],na.cc[3],na[10:12])

    output <- NULL
    output$B <- B
    output$beta <- beta
    output$X <- X
    output$N <- N
    output$strata <- strata
    output$n0 <- n0
    output$n1 <- n1
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
    if(cohort!=TRUE){
      output$NI <- NI
    }
    output$mean <- matrix.mean
    output$bias.mean <- matrix.bias
    output$pct.bias.mean <- matrix.pctbias
    output$median <- matrix.med
    output$bias.median <- matrix.bias.med
    output$pct.bias.med <- matrix.pctbias.med
    output$sd <- matrix.sd
    output$relative.uncertainty <- matrix.relunc
    output$mean.squared.error <- matrix.mse
    output$reported.standard.error <- matrix.sdave
    output$sd.reported.vs.actual <- matrix.rel.sd
    output$coverage.probability <- matrix.covprob
    output$power <- matrix.power
    output$na <- matrix.na

    class(output) <- "tpsSim"
    return(output)
     
  }
  else{ # CC
    if(length(n0)*length(n1)!=1){
      print("Error: 'n0' and 'n1' are not single numbers")
      return(-1)
    }
    cc.same <- 0
    if(n0==n1&is.null(ccDesign[1])==1){
      ##when balanced CC is specified, referent group should be specified
      cc.same <- 1
      ##print("Error: referent of the same balanced case-control designs")
      ##return(-1)
    }
    if(is.null(ccDesign[1])!=1){
      if(n0==ccDesign[1] && n1==ccDesign[2]){
        print("Warning: referent of the same case-control designs")
        cc.same <- 1
      }
    }
    if(!is.null(NI[1])){
      cohort <- FALSE
    }
    
    ##simulation 
    result <- mat.or.vec(nr=dim.beta*4,nc=B)
    for(i in 1:B){
      if(i %% monitor==0)
        cat("Repetition",i,"of",B,"complete\n") 
      result[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=1,n0=n0,n1=n1,cohort=cohort,NI=NI)
    }
    
    ##simulation for case-control
    result.cc <- mat.or.vec(nr=dim.beta*4,nc=B)
    if(is.null(ccDesign[1])==1){ #not specified by user, default is balanced CC
      ccDesign <- rep(times=2, round(sum(n0,n1)/2))
    }
    for(i in 1:B){
      if(i %% monitor==0)
        cat("Repetition",i,"of",B,"for a referent case-control design complete\n") 
      result.cc[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=1,n0=ccDesign[1],n1=ccDesign[2],cohort=cohort,NI=NI)
    }
    
    ##NA
    na.I <- numeric(3)
    na <- numeric(3)
    na.cc <- numeric(3)

    ##NA by too large deviation of estimates from truth in the second element of beta
    ## possibly due to numerical difficulty
    na.I[2] <- sum(apply(result[(2*dim.beta+1):(3*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
    na[2] <- sum(apply(result[2:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta[-1],threshold=threshold)) #intercept estimate is not reliable in case-control
    na.cc[2] <- sum(apply(result.cc[2:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta[-1],threshold=threshold)) #intercept estimate is not reliable in case-control   
 
    ##NA by sigularity of design matrix
    na.I[3] <- sum(is.na(result[(3*dim.beta+1),]))
    na[3] <- sum(is.na(result[(dim.beta+1),]))
    na.cc[3] <- sum(is.na(result.cc[(dim.beta+1),]))

    ##NA for all estimates in each method
    na.I[1] <- sum(is.na(result[(2*dim.beta+1),]))
    na[1] <- sum(is.na(result[1,]))
    na.cc[1] <- sum(is.na(result.cc[1,]))
    
    ##NA removed
    result.original <- result
    result.cc.original <- result.cc
    for(i in 1:dim.beta){
      result[c((2*dim.beta+1):(4*dim.beta)),
             result[(2*dim.beta+i),]-beta[i] < threshold[1]
             |result[(2*dim.beta+i),]-beta[i] > threshold[2]
             | is.na(result[(3*dim.beta+1),])
             | is.na(result[(2*dim.beta+1),])] <- NA
      if(i!=1){
        result[c(1:(2*dim.beta)),
               result[i,]-beta[i] < threshold[1]
               |result[i,]-beta[i] > threshold[2]
               | is.na(result[(dim.beta+1),])
               | is.na(result[1,])] <- NA
        result.cc[c(1:(2*dim.beta)),
                  result.cc[i,]-beta[i] < threshold[1]
                  |result.cc[i,]-beta[i] > threshold[2]
                  | is.na(result.cc[(dim.beta+1),])
                  | is.na(result.cc[1,])] <- NA
      }
    }

    ##mean
    est.ave.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,mean,na.rm=TRUE)
    est.ave <- apply(result[1:dim.beta,],1,mean,na.rm=TRUE)
    est.ave.cc <- apply(result.cc[1:dim.beta,],1,mean,na.rm=TRUE)

    ##(percent) bias
    bias.I <- est.ave.I-beta
    bias <- est.ave-beta
    bias.cc <- est.ave.cc-beta

    pct.bias.I <- (est.ave.I-beta)/abs(beta)*100
    pct.bias <- (est.ave-beta)/abs(beta)*100
    pct.bias.cc <- (est.ave.cc-beta)/abs(beta)*100
    beta.0.ind <- which(beta==0) #consideration of case of dividing by zero
    if(!length(beta.0.ind)){
      num.null.0 <- dim.beta
    }else{  
      num.null.0 <- dim.beta - length(beta.0.ind)
    }
    if(length(beta.0.ind) > 0){
      pct.bias.I <- pct.bias.I[-beta.0.ind]
      pct.bias <- pct.bias[-beta.0.ind]
      pct.bias.cc <- pct.bias.cc[-beta.0.ind]
    }

    ##median
    med.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,median,na.rm=TRUE)
    med <- apply(result[1:dim.beta,],1,median,na.rm=TRUE)
    med.cc <- apply(result.cc[1:dim.beta,],1,median,na.rm=TRUE)

    ##percent bias
    pct.bias.med.I <- (med.I-beta)/abs(beta)*100
    pct.bias.med <- (med-beta)/abs(beta)*100
    pct.bias.med.cc <- (med.cc-beta)/abs(beta)*100
    if(length(beta.0.ind) > 0){
      pct.bias.med.I <- pct.bias.med.I[-beta.0.ind]
      pct.bias.med <- pct.bias.med[-beta.0.ind]
      pct.bias.med.cc <- pct.bias.med.cc[-beta.0.ind]
    }
    
    ##standard deviation
    sds.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,sd,na.rm=TRUE)
    sds <- apply(result[1:dim.beta,],1,sd,na.rm=TRUE)
    sds.cc <- apply(result.cc[1:dim.beta,],1,sd,na.rm=TRUE)
     
    ##relative uncertainty
    if(referent==1){
      rel.unc.I <- rep(100.000,times=dim.beta)
      rel.unc <- sds/sds.I*100
      rel.unc.cc <- sds.cc/sds.I*100
    }else if(referent==2){
      rel.unc.I <- sds.I/sds.cc*100
      rel.unc <- sds/sds.cc*100
      rel.unc.cc <- rep(100.000,times=dim.beta)
    }else if(referent==3){
      rel.unc.I <- sds.I/sds*100
      if(cc.same==0){
        rel.unc <- rep(100.000,times=dim.beta)
        rel.unc.cc <- sds.cc/sds*100
      }else if(cc.same==1){
        rel.unc <- rep(100.000,times=dim.beta)
        rel.unc.cc <- rep(100.000,times=dim.beta)
      }
    }
      
    
    ##sd average
    sd.ave.I <- apply(result[(3*dim.beta+1):(4*dim.beta),],1,mean,na.rm=TRUE)
    sd.ave <- apply(result[(dim.beta+1):(2*dim.beta),],1,mean,na.rm=TRUE)
    sd.ave.cc <- apply(result.cc[(dim.beta+1):(2*dim.beta),],1,mean,na.rm=TRUE)

    ##ratio of reported sds and the actual truth
    rel.sd.I <- sd.ave.I/sds.I
    rel.sd <- sd.ave/sds
    rel.sd.cc <- sd.ave.cc/sds.cc
    
    ##MSE
    mse.I <- (est.ave.I-beta)^2+sds.I^2
    mse <- (est.ave-beta)^2+sds^2
    mse.cc <- (est.ave.cc-beta)^2+sds.cc^2
    
    ##coverage probability
    cov.prob.I <- 100*apply(ifelse(result[(2*dim.beta+1):(3*dim.beta),] + qnorm(alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                             < beta &
                             beta < result[(2*dim.beta+1):(3*dim.beta),]+qnorm(1-alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                             ,1,0),1,mean,na.rm=TRUE)

    cov.prob <- 100*apply(ifelse(result[1:dim.beta,] + qnorm(alpha/2)*result[(dim.beta+1):(2*dim.beta),]
                              < beta &
                           beta < result[1:dim.beta,]+qnorm(1-alpha/2)*result[(dim.beta+1):(2*dim.beta),]
                           ,1,0),1,mean,na.rm=TRUE)
  
    cov.prob.cc <- 100*apply(ifelse(result.cc[1:dim.beta,] + qnorm(alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                              < beta &
                           beta < result.cc[1:dim.beta,]+qnorm(1-alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                           ,1,0),1,mean,na.rm=TRUE)

    ##power with H0:beta=0
    power.0.I <- 100-100*apply(ifelse(result[(2*dim.beta+1):(3*dim.beta),]+ qnorm(alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                              < 0 &
                              0 < result[(2*dim.beta+1):(3*dim.beta),]+qnorm(1-alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                              ,1,0),1,mean,na.rm=TRUE)

    power.0 <- 100-100*apply(ifelse(result[1:dim.beta,]+ qnorm(alpha/2)*result[(dim.beta+1):(2*dim.beta),]
                                 < 0 &
                                 0 < result[1:dim.beta,]+qnorm(1-alpha/2)*result[(dim.beta+1):(2*dim.beta),]
                                 ,1,0),1,mean,na.rm=TRUE)

    power.0.cc <- 100-100*apply(ifelse(result.cc[1:dim.beta,]+ qnorm(alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                                 < 0 &
                                 0 < result.cc[1:dim.beta,]+qnorm(1-alpha/2)*result.cc[(dim.beta+1):(2*dim.beta),]
                                 ,1,0),1,mean,na.rm=TRUE)
    if(length(beta.0.ind) > 0){
      power.0.I <- power.0.I[-beta.0.ind]
      power.0 <- power.0[-beta.0.ind]
      power.0.cc <- power.0.cc[-beta.0.ind]
    }

    ## create matrices for output
    cov.name <- betaNames
    
    matrix.mean <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.mean) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.mean) <- cov.name
    matrix.mean[1,] <- est.ave.I
    matrix.mean[3,] <- est.ave
    matrix.mean[2,] <- est.ave.cc
    matrix.mean[c(2,3),1] <- NA
    
    matrix.bias <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.bias) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.bias) <- cov.name
    matrix.bias[1,] <- bias.I
    matrix.bias[3,] <- bias
    matrix.bias[2,] <- bias.cc
    matrix.bias[c(2,3),1] <- NA
    
    if(num.null.0!=0){
      matrix.pctbias <- as.matrix(mat.or.vec(nc=num.null.0,nr=3))
      rownames(matrix.pctbias) <- c("CD","C-C(ref)","Case-Control")
      if(!length(beta.0.ind)){
        colnames(matrix.pctbias) <- cov.name
      }else{
        colnames(matrix.pctbias) <- cov.name[-beta.0.ind]
      }
      matrix.pctbias[1,] <- pct.bias.I[1:num.null.0]
      matrix.pctbias[3,] <- pct.bias[1:num.null.0]  
      matrix.pctbias[2,] <- pct.bias.cc[1:num.null.0]  
      matrix.pctbias[c(2,3),1] <- NA
    }

    matrix.med <- mat.or.vec(nc=dim.beta,nr=3)
    matrix.med[1,] <- med.I
    matrix.med[3,] <- med
    matrix.med[2,] <- med.cc
    matrix.med[c(2,3),1] <- NA

    matrix.bias.med <- matrix.med
    matrix.beta <- rbind(beta,beta,beta)
    matrix.bias.med <- matrix.bias.med-matrix.beta
    
    if(num.null.0!=0){
      matrix.pctbias.med <- as.matrix(mat.or.vec(nc=num.null.0,nr=3))
      rownames(matrix.pctbias.med) <- c("CD","C-C(ref)","Case-Control")
      if(!length(beta.0.ind)){
        colnames(matrix.pctbias.med) <- cov.name
      }else{
        colnames(matrix.pctbias.med) <- cov.name[-beta.0.ind]
      }
      matrix.pctbias.med[1,] <- pct.bias.med.I[1:num.null.0]
      matrix.pctbias.med[3,] <- pct.bias.med[1:num.null.0]  
      matrix.pctbias.med[2,] <- pct.bias.med.cc[1:num.null.0]  
      matrix.pctbias.med[c(2,3),1] <- NA
    }
    
    matrix.sd <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.sd) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.sd) <- cov.name
    matrix.sd[1,] <- sds.I
    matrix.sd[3,] <- sds
    matrix.sd[2,] <- sds.cc
    matrix.sd[c(2,3),1] <- NA
    
    matrix.relunc <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.relunc) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.relunc) <- cov.name
    matrix.relunc[1,] <- rel.unc.I
    matrix.relunc[3,] <- rel.unc
    matrix.relunc[2,] <- rel.unc.cc
    matrix.relunc[c(2,3),1] <- NA
    
    matrix.mse <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.mse) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.mse) <- cov.name
    matrix.mse[1,] <- mse.I
    matrix.mse[3,] <- mse
    matrix.mse[2,] <- mse.cc
    matrix.mse[c(2,3),1] <- NA
    
    matrix.sdave <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.sdave) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.sdave) <- cov.name
    matrix.sdave[1,] <- sd.ave.I
    matrix.sdave[3,] <- sd.ave
    matrix.sdave[2,] <- sd.ave.cc
    matrix.sdave[c(2,3),1] <- NA
    
    matrix.rel.sd <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.rel.sd) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.rel.sd) <- cov.name
    matrix.rel.sd[1,] <- rel.sd.I
    matrix.rel.sd[3,] <- rel.sd
    matrix.rel.sd[2,] <- rel.sd.cc
    matrix.rel.sd[c(2,3),1] <- NA
    
    matrix.covprob <- mat.or.vec(nc=dim.beta,nr=3)
    rownames(matrix.covprob) <- c("CD","C-C(ref)","Case-Control")
    colnames(matrix.covprob) <- cov.name
    matrix.covprob[1,] <- cov.prob.I
    matrix.covprob[3,] <- cov.prob
    matrix.covprob[2,] <- cov.prob.cc
    matrix.covprob[c(2,3),1] <- NA
    
    if(num.null.0!=0){
      matrix.power <- as.matrix(mat.or.vec(nc=num.null.0,nr=3))
      rownames(matrix.power) <- c("CD","C-C(ref)","Case-Control")
      if(!length(beta.0.ind)){
        colnames(matrix.power) <- cov.name
      }else{
        colnames(matrix.power) <- cov.name[-beta.0.ind]
      }
      matrix.power[1,] <- power.0.I[1:num.null.0]
      matrix.power[3,] <- power.0[1:num.null.0]
      matrix.power[2,] <- power.0.cc[1:num.null.0]
      matrix.power[c(2,3),1] <- NA
    }
    
    matrix.na <- mat.or.vec(nc=3,nr=3)
    colnames(matrix.na) <- c("CD","C-C(ref)","Case-Control")
    rownames(matrix.na) <- c("Point Estimate","Above Threshold","Standard Deviation")
    matrix.na[1,] <- c(na.I[1],na[1],na.cc[1])
    matrix.na[3,] <- c(na.I[2],na[2],na.cc[2])
    matrix.na[2,] <- c(na.I[3],na[3],na.cc[3])

    if(cc.same==1){
      matrix.mean <- matrix.mean[c(1,3),]
      matrix.bias <- matrix.bias[c(1,3),]
      matrix.pctbias <- matrix.pctbias[c(1,3),]
      matrix.med <- matrix.med[c(1,3),]
      matrix.bias.med <- matrix.bias.med[c(1,3),]
      matrix.pctbias.med <- matrix.pctbias.med[c(1,3),]
      matrix.sd <- matrix.sd[c(1,3),]
      matrix.relunc <- matrix.relunc[c(1,3),]
      matrix.mse <- matrix.mse[c(1,3),]
      matrix.sdave <- matrix.sdave[c(1,3),]
      matrix.rel.sd <- matrix.rel.sd[c(1,3),]
      matrix.covprob <- matrix.covprob[c(1,3),]
      matrix.power <- matrix.power[c(1,3),]
      matrix.na <- matrix.na[,c(1,3)]
    }
    
    output <- NULL
    output$B <- B
    output$beta <- beta
    output$X <- X
    output$N <- N
    output$strata <- strata
    output$n0 <- n0
    output$n1 <- n1
    output$ccDesign <- ccDesign
    output$referent <- referent
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
    output$pct.bias.mean <- matrix.pctbias
    output$median <- matrix.med
    output$bias.median <- matrix.bias.med
    output$pct.bias.med <- matrix.pctbias.med
    output$sd <- matrix.sd
    output$relative.uncertainty <- matrix.relunc 
    output$mean.squared.error <- matrix.mse
    output$reported.standard.error <- matrix.sdave
    output$sd.reported.vs.actual <- matrix.rel.sd
    output$coverage.probability <- matrix.covprob
    output$power <- matrix.power
    output$na <- matrix.na

    class(output) <- "tpsSim"
    
    return(output)
  }
}

