rmvhyper <-
function(Mk, m)
{
  K <- length(Mk)
  M <- sum(Mk)
  mk <- rep(0, K)
	
  ##	
  if(m > M){
    mk <- Mk
  }else if(K==1){
    mk <- m
  }else{
    mk[1] <- rhyper(1, Mk[1], M - Mk[1], m)
    if(K > 2){
      for(j in 2:(K-1)){
        mk[j] <- rhyper(1, Mk[j], M - sum(Mk[1:j]), m - sum(mk[1:(j-1)]))
      }
    }
    mk[K] <- m - sum(mk[1:(K-1)])
  }
  return(mk)
}

