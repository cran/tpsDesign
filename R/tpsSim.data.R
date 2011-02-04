tpsSim.data <-
function(betaTruth, X, N, strata, n0, n1, cohort=TRUE, NI=NULL){

  beta <- betaTruth
  strata <- sort(strata)
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
  Death <- rbinom(n=rep(x=1,times=length(N)),size=N, prob=prob)
  Alive <- N-Death

  ## case-control option
  if(cohort!=TRUE){
    Alive <- round(Alive/sum(Alive)*NI[1])
    Death <- round(Death/sum(Death)*NI[2]) 
  }
  
  ## stratified indicator
  ## typically 1*level at X1 + 10*level at X2 etc.
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
      
  ## match up strata indicator and sample size assignment using reverse lexicographic order
  ## order of (0 10 20 1 11 21) is ( 1 2 3 4 5 6)
  lev <- unique(sort(st.ind))
  lev <- lev[order(strReverse(as.character(lev)))] 
      
  ## Phase 1 data
  N.cont <- tapply(Alive, st.ind, FUN=sum)
  N.case <- tapply(Death, st.ind, FUN=sum)
  
  ## Phase 2 data
  AliveII <- Alive
  DeathII <- Death
  for(i in 1:length(unique(st.ind))){
    ##
    AliveII[st.ind == lev[i]] <- rmvhyper(Alive[st.ind == lev[i]], n0[i])
    DeathII[st.ind == lev[i]] <- rmvhyper(Death[st.ind == lev[i]], n1[i])
  }
  
  return(list(AliveI = Alive, DeathI=Death, AliveII = as.integer(AliveII),
              DeathII = as.integer(DeathII), N.control = N.cont, N.case = N.case,
              strata.index=st.ind))
}

