coreChecks <-
function(betaTruth,
												X,
												N,
												betaNames)
{
	##
	if(ncol(as.matrix(X)) == 1)
    return("* 'X' should have at least two columns")
	
	##
	if(sum(X[,1] != 1) > 0)
  	return("* 'X' requires the first column to represent the intercept")
  
	##
	for(i in 2:ncol(X))
	{
		prob1 <- sum(ceiling(X[,i]) != floor(X[,i]))
		prob2 <- (min(X[,i]) != 0)
		prob3 <- (max(X[,i]) != (length(unique(X[,i])) - 1))
		prob4 <- (max(X[,i]) == 0)
		if((prob1 + prob2 + prob3 + prob4) > 0)
			return("* check that each variable is consistent with the {0,1,2,...} coding convention")
	}
	
  ##
  if(length(N) != nrow(X))
    return("* incompatible dimensions of 'X' and 'N'")
	
  ##
  p <- sum(unlist(lapply(apply(X, 2, unique), FUN=length)) - 1) + 1
  if(length(betaTruth) != p)
    return("* invalid dimension of 'beta'")

	##
	if(!is.null(betaNames))
	{
		if(length(betaTruth) != length(betaNames))
			return("* 'beta' and 'betaNames' are not of the same length")
		if(!is.character(betaNames))
		return("* elements of 'betaNames' are not character")
	}
  
	##
	return("")
}

