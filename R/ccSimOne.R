ccSimOne <-
function(B=1000, betaTruth, X, N, n0, n1,
                    alpha=.05, threshold=c(-Inf,Inf),monitor=NULL,
                    NI=NULL){
  
  beta <- betaTruth
  dim.beta <- length(beta)
    
  if(!is.null(NI[1])){
    cohort <- FALSE
  }else{
    cohort <- TRUE
  }
  
  ##simulation 
  result <- mat.or.vec(nr=dim.beta*4,nc=B)
  for(i in 1:B){
    if(i %% monitor==0)
      cat("Repetition",i,"of",B,"complete\n") 
    result[,i] <- tpsSim.fit(betaTruth=beta,X=X,N=N,strata=1,n0=n0,n1=n1,cohort=cohort,NI=NI)
  }
  
  
  ##NA
  na.I <- numeric(3)
  na <- numeric(3)
  
  ##NA by too large deviation of estimates from truth in the second element of beta
  ## possibly due to numerical difficulty
  na.I[2] <- sum(apply(result[(2*dim.beta+1):(3*dim.beta),], 2, FUN=tpsSim.na.big, betaTruth=beta,threshold=threshold))
  na[2] <- sum(apply(result[2:dim.beta,], 2, FUN=tpsSim.na.big, betaTruth=beta[-1],threshold=threshold)) #intercept estimate is not reliable in case-control
  
  ##NA by sigularity of design matrix
  na.I[3] <- sum(is.na(result[(3*dim.beta+1),]))
  na[3] <- sum(is.na(result[(dim.beta+1),]))

  ##NA for all estimates in each method
  na.I[1] <- sum(is.na(result[(2*dim.beta+1),]))
  na[1] <- sum(is.na(result[1,]))
    
  ##NA removed
  result.original <- result
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
    }
  }
  
  ##mean
  est.ave.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,mean,na.rm=TRUE)
  est.ave <- apply(result[1:dim.beta,],1,mean,na.rm=TRUE)
 
  ##(percent) bias
  bias.I <- est.ave.I-beta
  bias <- est.ave-beta
   
  pct.bias.I <- (est.ave.I-beta)/abs(beta)*100
  pct.bias <- (est.ave-beta)/abs(beta)*100
   beta.0.ind <- which(beta==0) #consideration of case of dividing by zero
  if(!length(beta.0.ind)){
    num.null.0 <- dim.beta
  }else{  
    num.null.0 <- dim.beta - length(beta.0.ind)
  }
  if(length(beta.0.ind) > 0){
    pct.bias.I <- pct.bias.I[-beta.0.ind]
    pct.bias <- pct.bias[-beta.0.ind]
  }
  
  ##median
  med.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,median,na.rm=TRUE)
  med <- apply(result[1:dim.beta,],1,median,na.rm=TRUE)
   
  ##percent bias
  pct.bias.med.I <- (med.I-beta)/abs(beta)*100
  pct.bias.med <- (med-beta)/abs(beta)*100
  if(length(beta.0.ind) > 0){
    pct.bias.med.I <- pct.bias.med.I[-beta.0.ind]
    pct.bias.med <- pct.bias.med[-beta.0.ind]
  }
  
  ##standard deviation
  sds.I <- apply(result[(2*dim.beta+1):(3*dim.beta),],1,sd,na.rm=TRUE)
  sds <- apply(result[1:dim.beta,],1,sd,na.rm=TRUE)
  
    
  ##sd average
  sd.ave.I <- apply(result[(3*dim.beta+1):(4*dim.beta),],1,mean,na.rm=TRUE)
  sd.ave <- apply(result[(dim.beta+1):(2*dim.beta),],1,mean,na.rm=TRUE)
  
  ##ratio of reported sds and the actual truth
  rel.sd.I <- sd.ave.I/sds.I
  rel.sd <- sd.ave/sds
   
  ##MSE
  mse.I <- (est.ave.I-beta)^2+sds.I^2
  mse <- (est.ave-beta)^2+sds^2
   
  ##coverage probability
  cov.prob.I <- 100*apply(ifelse(result[(2*dim.beta+1):(3*dim.beta),] + qnorm(alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                                 < beta &
                                 beta < result[(2*dim.beta+1):(3*dim.beta),]+qnorm(1-alpha/2)*result[(3*dim.beta+1):(4*dim.beta),]
                                 ,1,0),1,mean,na.rm=TRUE)
  
  cov.prob <- 100*apply(ifelse(result[1:dim.beta,] + qnorm(alpha/2)*result[(dim.beta+1):(2*dim.beta),]
                               < beta &
                               beta < result[1:dim.beta,]+qnorm(1-alpha/2)*result[(dim.beta+1):(2*dim.beta),]
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
  
  if(length(beta.0.ind) > 0){
    power.0.I <- power.0.I[-beta.0.ind]
    power.0 <- power.0[-beta.0.ind]
  }

  ## create matrices for output
  
  matrix.mean <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.mean[1,] <- est.ave.I
  matrix.mean[2,] <- est.ave
  matrix.mean[2,1] <- NA
  
  matrix.bias <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.bias[1,] <- bias.I
  matrix.bias[2,] <- bias
  matrix.bias[2,1] <- NA
  
  if(num.null.0!=0){
    matrix.pctbias <- mat.or.vec(nc=num.null.0,nr=2)
    matrix.pctbias[1,] <- pct.bias.I[1:num.null.0]
    matrix.pctbias[2,] <- pct.bias[1:num.null.0]  
    matrix.pctbias[2,1] <- NA
  }
  
  matrix.med <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.med[1,] <- med.I
  matrix.med[2,] <- med
  matrix.med[2,1] <- NA
  
  matrix.bias.med <- matrix.med
  matrix.beta <- rbind(beta,beta)
  matrix.bias.med <- matrix.bias.med-matrix.beta
    
  if(num.null.0!=0){
    matrix.pctbias.med <- mat.or.vec(nc=num.null.0,nr=2)
    matrix.pctbias.med[1,] <- pct.bias.med.I[1:num.null.0]
    matrix.pctbias.med[2,] <- pct.bias.med[1:num.null.0]  
    matrix.pctbias.med[2,1] <- NA
  }
    
  matrix.sd <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.sd[1,] <- sds.I
  matrix.sd[2,] <- sds
  matrix.sd[2,1] <- NA
    
  matrix.mse <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.mse[1,] <- mse.I
  matrix.mse[2,] <- mse
  matrix.mse[2,1] <- NA
  
  matrix.sdave <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.sdave[1,] <- sd.ave.I
  matrix.sdave[2,] <- sd.ave
  matrix.sdave[2,1] <- NA
  
  matrix.rel.sd <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.rel.sd[1,] <- rel.sd.I
  matrix.rel.sd[2,] <- rel.sd
  matrix.rel.sd[2,1] <- NA
  
  matrix.covprob <- mat.or.vec(nc=dim.beta,nr=2)
  matrix.covprob[1,] <- cov.prob.I
  matrix.covprob[2,] <- cov.prob
  matrix.covprob[2,1] <- NA
  
  if(num.null.0!=0){
    matrix.power <- mat.or.vec(nc=num.null.0,nr=2)
    matrix.power[1,] <- power.0.I[1:num.null.0]
    matrix.power[2,] <- power.0[1:num.null.0]
    matrix.power[2,1] <- NA
  }
  
  matrix.na <- mat.or.vec(nc=2,nr=3)
  matrix.na[1,] <- c(na.I[1],na[1])
  matrix.na[2,] <- c(na.I[2],na[2])
  matrix.na[3,] <- c(na.I[3],na[3])
  
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

