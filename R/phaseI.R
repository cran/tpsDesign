phaseI <-
function(betaTruth, X, N, strata=NULL, n0=NULL, n1=NULL,digits=NULL){
  
  beta <- betaTruth

  ##check compatibility
  if(length(N)!=nrow(X)){
    print("Error: invalid dimensions of 'N' and 'X'")
    return(-1)
  }
  
  ## #columns of design matrix
  n.col <- ncol(X)
  
  ## #levels in each covariates
  n.lev <- numeric((n.col-1)) # intercept excluded
  for(i in 1 : (n.col-1) ){
    n.lev[i] <- length(unique(X[ , (i+1)] ))
  }
  n.lev.cum <- c(0,cumsum(n.lev-1)) #in all 1st groups =baseline(not counted)
  if(length(beta)!=rev(n.lev.cum+1)[1]){
    print("Error: invalid dimension of 'beta'")
    return(-1)
  }
  
  ##matrix of coefficients
  P <- mat.or.vec(nc=n.col,nr=nrow(X)) #add intercept column later
  for(i in 2:n.col){   # i : column
    for(j in 1:(n.lev[(i-1)]-1)){  # j : indicator of level
      P[X[,i]==j,i] <- beta[(1+n.lev.cum[(i-1)]+j)]
    }
  }
  P[,1] <- beta[1]
  
  ##generating data
  prob <- exp(rowSums(P))/(1+exp(rowSums(P)))
  exp.death <- round(N*prob)
  exp.alive <- N-exp.death

  ## only phaseI information when either strata, n0, or n1 is null
  if(is.null(strata)){
    ## naming the row names of matrices
    nrows <- nrow(X)
    phaseI <- mat.or.vec(nc=2,nr=nrows)
    colnames(phaseI) <- c("Y=0","Y=1")
    rownames(phaseI) <- rep(times=nrows,"a")
    phaseI[,1] <- exp.alive
    phaseI[,2] <- exp.death
    
    for(i in 1:nrows){
      rownames(phaseI)[i] <- as.character(paste(colnames(X[,2:ncol(X)]),"=",X[i,2:ncol(X)]," ",collapse=" "))
    }
    cat("Expected Phase I:\n")
    print(phaseI,digits=digits)
    cat("\n")
    
    if(!is.null(n0)){
      print("Warning: 'strata' is NULL")
      if(is.null(n1)){
        print("Warning: 'n1' is NULL")
      }
    }else if(!is.null(n1)){
      print("Warning: 'strata' is NULL")
           print("Warning: 'n0' is NULL")
    }
    invisible()
  }else{
    ## check strata
    if(length(intersect(strata,1:ncol(X)))==0){
      print("Error: 'strata' is invalid")
      invisible()
      return()
    }
    strata <- sort(strata)
    
    ##collapse strata
    st.ind <- X[,strata[1]]
    if(strata[1]!=1){
      if(length(strata)>1){
        j <- 1 #power for base 10
        if(n.lev[(strata[1]-1)]>9){
          st.ind <- st.ind*10
          j <- 2
        }
        for(i in 2:length(strata)){
          if(n.lev[(strata[i]-1)]> 9) j <- j+1 
          st.ind <- st.ind + 10^j * X[,strata[i]]
          j <- j+1
      }
      }
    }  
    temp <- st.ind
    lev <- unique(sort(st.ind))
    lev <- lev[order(strReverse(as.character(lev)))] 
    for(i in 1:length(lev)) st.ind[temp == lev[i]] <- i
    
    
    ## naming the row names of matrices
    nrows <- prod(n.lev[(strata-1)])
    phaseI <- mat.or.vec(nc=2,nr=nrows)
    colnames(phaseI) <- c("Y=0","Y=1")
    rownames(phaseI) <- rep(times=nrows,"a")
    if(length(strata)>1){ ##stratified by more than one stratum   
      mat.st.ind <- mat.or.vec(nr=nrows,nc=length(strata))
      colnames(mat.st.ind) <- colnames(X)[strata]
      lev.length <- n.lev[(strata[1]-1)]
      n.rep <- nrows/lev.length
      
      n.lev.st <- n.lev[(strata-1)]
      for(i in 1:length(strata)){
        s <- NULL
        for(j in 1:n.lev[(strata[i]-1)]){
        s <- c(s,rep(times=n.rep,j-1))
      }
        lev.length <- length(s)
        if(i!=length(strata)){
          n.rep <- n.rep/n.lev.st[(i+1)]
        }
      mat.st.ind[,i] <- rep(times=nrows/lev.length,s)
      }
      
      for(i in 1:nrows){
        rownames(phaseI)[i] <- as.character(paste(colnames(mat.st.ind),"=",mat.st.ind[i,]," ",collapse=" "))
      }
    }else{ ##stratified by one stratum
      for(i in 1:nrows){
        rownames(phaseI)[i] <- as.character(paste(colnames(X)[strata],"=",i))
      }
    }
    
    for(i in 1:nrows){
      phaseI[i,1] <- sum(exp.alive[st.ind==i])
      phaseI[i,2] <- sum(exp.death[st.ind==i])
    }
    
    cat("Expected Phase I:\n")
    print(phaseI,digits=digits)
    cat("\n")

    if(!is.null(n0)&&!is.null(n1)){

      ##check dimesion of n1 and n0
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

      phaseII <- phaseI
      phaseII[,1] <- n0
      phaseII[,2] <- n1
      
      a <- min(phaseII/phaseI)
      for(i in 1:100){
        if(10^(-i)<a)
          break
      }
      if(is.null(digits)){
        digits <- i+2
      }
      
      phaseIIProb <- round(phaseII/phaseI, digits=digits)
      
      
      cat("Phase II:\n")
      print(phaseII,digits=digits)
      cat("\n")
      
      cat("Expected Phase II sampling probabilities:\n")
      print(round(phaseIIProb,digits=digits))
      cat("\n")

      for(i in 1:nrow(phaseI)){
        if(phaseI[i,1]<phaseII[i,1]){
          cat("Warning: strata", rownames(phaseI)[i],"does not have enough controls\n")
        }
        if(phaseI[i,2]<phaseII[i,2]){
          cat("Warning: strata", rownames(phaseI)[i],"does not have enough cases\n")
        }
      }
      invisible()
      
    }
    if(is.null(n0)){
      print("Warning: 'n0' is NULL")
    }
    if(is.null(n1)){
      print("Warning: 'n1' is NULL")
    }
  }
  invisible()
}

