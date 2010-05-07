tpsSim.fit <-
function(betaTruth, X, N, strata, n0, n1, var.name=NULL){

  beta <- betaTruth
  strata <- sort(strata)
  ##simulate data
  data.sim <- tpsSim.data(betaTruth=beta,X=X,N=N,strata=strata,n0=n0,n1=n1)
  if(data.sim[[1]][1]==-1) return(-1)
  N.control <- data.sim$N.control
  N.case <- data.sim$N.case
  
  ##formula and rename design matrix
  p <- ncol(X)
  if(is.null(var.name)){
    formI <- "cbind(Death,Alive) ~ factor(X1)"
    formII <- "cbind(DeathII,AliveII) ~ factor(X1)"
    if(p > 2){
      for(i in 2: (p-1)){
        formI <- paste(formI, "+ factor(X", i, ")", sep="")
        formII <- paste(formII, "+ factor(X", i, ")", sep="")
      }
    }
    colnames(X) <- c("Intercept",paste("X",1:(p-1),sep=""))
  }else{
    formI <- paste("cbind(Death,Alive) ~ factor(",var.name[1],")",sep="")
    
    formII <- paste("cbind(DeathII,AliveII) ~ factor(",var.name[1],")",sep="")
    if(p > 2){ ## at least 2 covariates
      for(i in 2: (p-1)){
        formI <- paste(formI, "+ factor(", var.name[i],")", sep="")
        formII <- paste(formII, "+ factor(", var.name[i],")", sep="")
      }
    }
  }

  
  ##create data.frame
  data.sim <- cbind(N, data.sim$strata.index, data.sim$DeathI, data.sim$AliveI,
                    data.sim$DeathII, data.sim$AliveII)
  colnames(data.sim) <- c("N", "Strata", "Death", "Alive", "DeathII", "AliveII")
  data.sim <- as.data.frame(cbind(X,data.sim))

  ## modify strata index for tps
  temp <- data.sim$Strata
  lev <- unique(sort(data.sim$Strata))
  for(i in 1:length(lev)) data.sim$Strata[temp == lev[i]] <- i

  if(length(lev)!=1){ #non CC
    ## estimates using all phase I data
    dim.beta <- length(beta)
    fit <- glm(formula=formI, data=data.sim, family=binomial)
    if(class(fit)[1]== "try-error"){
      beta.hat <- rep(NA, dim.beta)
      beta.hat.sd <- rep(NA, dim.beta)
    }
    if(class(fit)[1] == "glm"){
      beta.hat <- summary(fit)$coefficient[,1]
      beta.hat.sd <- summary(fit)$coefficient[,2]
    }
    ## WL
    op <- options()
    options(warn=-1)
    fitWL <- try(tps(formula=formII, data=data.sim, nn0=N.control, nn1=N.case,
                     group=data.sim$Strata, method="WL"),silent=TRUE)
    options(op)
    if(class(fitWL) == "try-error"){
      beta.hat.WL <- rep(NA, dim.beta)
      beta.hat.sd.WL <- rep(NA, dim.beta)
    }
    if(class(fitWL) == "tps"){
      beta.hat.WL <- fitWL$coef
      beta.hat.sd.WL <- sqrt(diag(fitWL$cove))
    }
    ## PL
    fitPL <- try(tps(formula=formII, data=data.sim, nn0=N.control, nn1=N.case,
                     group=data.sim$Strata, method="PL"),silent=TRUE)  
    if(class(fitPL) == "try-error"){
      beta.hat.PL <- rep(NA, dim.beta)
      beta.hat.sd.PL <- rep(NA, dim.beta)
    } 
    if(class(fitPL) == "tps"){
      beta.hat.PL <- fitPL$coef
      beta.hat.sd.PL <- sqrt(diag(fitPL$covm))
    }
    ## ML
    fitML <- try(tps(formula=formII, data=data.sim, nn0=N.control, nn1=N.case,
                     group=data.sim$Strata, method="ML"),silent=TRUE)  
    if(class(fitML) == "try-error"){
      beta.hat.ML <- rep(NA, dim.beta)
      beta.hat.sd.ML <- rep(NA, dim.beta)
    } 
    if(class(fitML) == "tps"){
      beta.hat.ML <- fitML$coef
      beta.hat.sd.ML <- sqrt(diag(fitML$covm))
    }
    
    return(c(beta.hat, beta.hat.WL, beta.hat.PL, beta.hat.ML,
             beta.hat.sd, beta.hat.sd.WL, beta.hat.sd.PL, beta.hat.sd.ML))
  }else{ # CC
    ## estimates using all phase I data
    dim.beta <- length(beta)
    fit <- glm(formula=formI, data=data.sim, family=binomial)
    if(class(fit)[1]== "try-error"){
      beta.hat <- rep(NA, dim.beta)
      beta.hat.sd <- rep(NA, dim.beta)
    }
    if(class(fit)[1] == "glm"){
      beta.hat <- summary(fit)$coefficient[,1]
      beta.hat.sd <- summary(fit)$coefficient[,2]
    }
    ## estimates for CC
    fitGLM <- try(glm(formula=formII, data=data.sim, family=binomial))
    if(class(fitGLM)[1] == "try-error"){
      beta.hat.cc <- rep(NA, dim.beta)
      beta.hat.sd.cc <- rep(NA, dim.beta)
    } 
    if(class(fitGLM)[1] == "glm"){
      beta.hat.cc <- fitGLM$coefficients
      beta.hat.sd.cc <- sqrt(diag(vcov(fitGLM)))
      if(length(beta.hat.cc)>length(beta.hat.sd.cc)){
        a <- is.na(beta.hat.cc)
        id <- which(a==TRUE)
        if(sum(id)==0)break
        beta.hat.sd.cc2 <- numeric(length(beta.hat.cc))
        beta.hat.sd.cc2[id] <- NA
        beta.hat.sd.cc2[-id] <- beta.hat.sd.cc
        beta.hat.sd.cc <- beta.hat.sd.cc2
      }
    }
    return(as.vector(c(beta.hat.cc, beta.hat.sd.cc, beta.hat, beta.hat.sd)))
  }
}

