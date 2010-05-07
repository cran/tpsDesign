tpsSimAllOne <-
function(B, betaTruth, X, N, strata, n0, n1,
                         alpha, threshold=c(-Inf,Inf), monitor=NULL){

  beta <- betaTruth
  ## valid strata
  if(length(intersect(strata,0:ncol(X)))==0){
    print("Error: 'strata' is invalid")
    return(-1)
  }

 
  dim.beta <- length(beta)
  strata <- sort(strata)
    
  ##simulation 
  result <- mat.or.vec(nr=dim.beta*8,nc=B)
  for(i in 1:B){
    if(i %% monitor==0)
      cat("Repetition",i,"of",B,"for a two-phase design complete\n") 
    result[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1)
  }
  ##NA
  na <- numeric(12)
      
  ##NA by too large deviation of estimates from truth possibly due to numerical difficulty
  na[5] <- sum(apply(result[1:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
  na[6] <- sum(apply(result[(dim.beta+1):(2*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
  na[7] <- sum(apply(result[(2*dim.beta+1):(3*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
  na[8] <- sum(apply(result[(3*dim.beta+1):(4*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
  
  
  ##NA by sigularity of design matrix
  na[9] <- sum(is.na(result[(4*dim.beta+1),]))
  na[10] <- sum(is.na(result[(5*dim.beta+1),]))
  na[11] <- sum(is.na(result[(6*dim.beta+1),]))
  for(i in 1:dim.beta){
    na[12] <- max(na[12],sum(is.na(result[(7*dim.beta+i),])))
  }
  
  ##NA for all estimates in each method
  na[1] <- sum(is.na(result[1,]))
  na[2] <- sum(is.na(result[(dim.beta+1),]))
  na[3] <- sum(is.na(result[(2*dim.beta+1),]))
  na[4] <- sum(is.na(result[(3*dim.beta+1),]))
  
  ##NA removed
  result.original <- result
  
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
    
  }
  
  ##mean
  est.ave <- apply(result[1:(4*dim.beta),],1,mean,na.rm=TRUE)
    
  ##percent bias in mean
  matrix.beta <- rbind(beta,beta,beta,beta)
  matrix.bias.mean <- matrix(est.ave,nrow=4,byrow=TRUE)-matrix.beta
  pct.bias <- matrix.bias.mean/abs(matrix.beta)*100
  beta.0.ind <- which(beta==0) #consideration of case of dividing by zero
  if(!length(beta.0.ind)){
    num.null.0 <- dim.beta
  }else{  
    num.null.0 <- dim.beta - length(beta.0.ind)
  }
  if(length(beta.0.ind) > 0){
    pct.bias <- pct.bias[,-beta.0.ind]
  }
 
  ##median
  med <- apply(result[1:(4*dim.beta),],1,median,na.rm=TRUE)
  

  ##percent bias in median
  matrix.diff <- matrix(med,nrow=4,byrow=TRUE)-matrix.beta
  pct.bias.med <- matrix.diff/abs(matrix.beta)*100
  if(length(beta.0.ind) > 0){
    pct.bias.med <- pct.bias.med[,-beta.0.ind]
  }    
    
  ##standard deviation
  sds <- apply(result[1:(4*dim.beta),],1,sd,na.rm=TRUE)
      
  ##sd average
  sd.ave <- apply(result[(4*dim.beta+1):(8*dim.beta),],1,mean,na.rm=TRUE)
  
  ##ratio of reported sds and the actual truth
  rel.sd <- sd.ave/sds
    
  ##MSE
  mse <- (est.ave-beta)^2+sds^2
      
  ##coverage probability
  cov.prob <- 100*apply(ifelse(result[1:(4*dim.beta),] + qnorm(alpha/2)*result[(4*dim.beta+1):(8*dim.beta),] < beta &
                               beta < result[1:(4*dim.beta),]+qnorm(1-alpha/2)*result[(4*dim.beta+1):(8*dim.beta),]
                               ,1,0),1,mean,na.rm=TRUE)
  
  ##power with H0:beta=0
  power.0 <- 100-100*apply(ifelse(result[1:(4*dim.beta),] + qnorm(alpha/2)*result[(4*dim.beta+1):(8*dim.beta),] < 0 &
                                  0 < result[1:(4*dim.beta),] + qnorm(1-alpha/2)*result[(4*dim.beta+1):(8*dim.beta),]
                                  ,1,0),1,mean,na.rm=TRUE)
  if(length(beta.0.ind) > 0){
    power.0 <- power.0[-c(beta.0.ind,(beta.0.ind+dim.beta),(beta.0.ind+2*dim.beta),(beta.0.ind+3*dim.beta))]
  }
    

      
  matrix.mean <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.mean[1,] <- est.ave[1:dim.beta]
  matrix.mean[2,] <- est.ave[(dim.beta+1):(2*dim.beta)]
  matrix.mean[3,] <- est.ave[(2*dim.beta+1):(3*dim.beta)]
  matrix.mean[4,] <- est.ave[(3*dim.beta+1):(4*dim.beta)]    
    
  matrix.bias <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.bias[1,] <- matrix.bias.mean[1,]
  matrix.bias[2,] <- matrix.bias.mean[2,]
  matrix.bias[3,] <- matrix.bias.mean[3,]
  matrix.bias[4,] <- matrix.bias.mean[4,]
  
  if(num.null.0!=0){
    matrix.pctbias <- mat.or.vec(nc=num.null.0,nr=4)
    matrix.pctbias[1,] <- pct.bias[1,]
    matrix.pctbias[2,] <- pct.bias[2,]
    matrix.pctbias[3,] <- pct.bias[3,]
    matrix.pctbias[4,] <- pct.bias[4,]
  }

  matrix.med <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.med[1,] <- med[1:dim.beta]
  matrix.med[2,] <- med[(dim.beta+1):(2*dim.beta)]
  matrix.med[3,] <- med[(2*dim.beta+1):(3*dim.beta)]
  matrix.med[4,] <- med[(3*dim.beta+1):(4*dim.beta)]    
  
  matrix.bias.med <- matrix.med
  matrix.beta <- rbind(beta,beta,beta,beta)
  matrix.bias.med <- matrix.bias.med-matrix.beta
    
  if(num.null.0!=0){
    matrix.pctbias.med <- mat.or.vec(nc=num.null.0,nr=4)
    matrix.pctbias.med[1,] <- pct.bias.med[1,]
    matrix.pctbias.med[2,] <- pct.bias.med[2,]
    matrix.pctbias.med[3,] <- pct.bias.med[3,]
    matrix.pctbias.med[4,] <- pct.bias.med[4,]
  }
  
  matrix.sd <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.sd[1,] <- sds[1:dim.beta]
  matrix.sd[2,] <- sds[(dim.beta+1):(2*dim.beta)]
  matrix.sd[3,] <- sds[(2*dim.beta+1):(3*dim.beta)]
  matrix.sd[4,] <- sds[(3*dim.beta+1):(4*dim.beta)]
      
  matrix.mse <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.mse[1,] <- mse[1:dim.beta]
  matrix.mse[2,] <- mse[(dim.beta+1):(2*dim.beta)]
  matrix.mse[3,] <- mse[(2*dim.beta+1):(3*dim.beta)]
  matrix.mse[4,] <- mse[(3*dim.beta+1):(4*dim.beta)]
      
  matrix.sdave <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.sdave[1,] <- sd.ave[1:dim.beta]
  matrix.sdave[2,] <- sd.ave[(dim.beta+1):(2*dim.beta)]
  matrix.sdave[3,] <- sd.ave[(2*dim.beta+1):(3*dim.beta)]
  matrix.sdave[4,] <- sd.ave[(3*dim.beta+1):(4*dim.beta)]
      
  matrix.rel.sd <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.rel.sd[1,] <- rel.sd[1:dim.beta]
  matrix.rel.sd[2,] <- rel.sd[(dim.beta+1):(2*dim.beta)]
  matrix.rel.sd[3,] <- rel.sd[(2*dim.beta+1):(3*dim.beta)]
  matrix.rel.sd[4,] <- rel.sd[(3*dim.beta+1):(4*dim.beta)]
      
  matrix.covprob <- mat.or.vec(nc=dim.beta,nr=4)
  matrix.covprob[1,] <- cov.prob[1:dim.beta]
  matrix.covprob[2,] <- cov.prob[(dim.beta+1):(2*dim.beta)]
  matrix.covprob[3,] <- cov.prob[(2*dim.beta+1):(3*dim.beta)]
  matrix.covprob[4,] <- cov.prob[(3*dim.beta+1):(4*dim.beta)]
      
  if(num.null.0!=0){
    matrix.power <- mat.or.vec(nc=num.null.0,nr=4)
    matrix.power[1,] <- power.0[1:num.null.0]
    matrix.power[2,] <- power.0[(num.null.0+1):(2*num.null.0)]
    matrix.power[3,] <- power.0[(2*num.null.0+1):(3*num.null.0)]
    matrix.power[4,] <- power.0[(3*num.null.0+1):(4*num.null.0)]
  }
    
  matrix.na <- mat.or.vec(nc=4,nr=3)
  matrix.na[1,] <- c(na[1],na[2:4])
  matrix.na[2,] <- c(na[5],na[6:8])
  matrix.na[3,] <- c(na[9],na[10:12])

  output <- NULL
  output$mean <- matrix.mean
  output$bias.mean <- matrix.bias
  output$pct.bias.mean <- matrix.pctbias
  output$median <- matrix.med
  output$bias.median <- matrix.bias.med
  output$pct.bias.med <- matrix.pctbias.med
  output$sd <- matrix.sd
  output$mean.squared.error <- matrix.mse
  output$reported.standard.error <- matrix.sdave
  output$sd.reported.vs.actual <- matrix.rel.sd
  output$coverage.probability <- matrix.covprob
  output$power <- matrix.power
  output$na <- matrix.na
  
  return(output)
}

