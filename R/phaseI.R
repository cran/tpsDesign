phaseI <-
function(betaTruth,
									 X,
									 N,
									 strata=NULL,
									 nII0=NULL,
									 nII1=NULL,
									 cohort=TRUE,
									 NI=NULL,
									 digits=NULL)
{
	## Checks
	##
	problem <- coreChecks(betaTruth=betaTruth, X=X, N=N, betaNames=NULL)
	if(problem != "")
		stop(problem)
 	##
	problem <- phaseIChecks(X=X, strata=strata, cohort=cohort, NI=NI, nII0=nII0, nII1=nII1)
	if(problem != "")
		stop(problem)

	##
  if(is.null(colnames(X)))
  	colnames(X) <- c("Int", paste("V", 1:(ncol(X) - 1), sep = ""))
	
  ## Calculate the expected number of Cases/nonCases for each row of X
  ##  * adjustment if the phase I sample is a case-control sample (of size NI=c(N0,N1))
  ##
  Xexp <- expandCatX(X)
  piY  <- expit(as.vector(Xexp %*% betaTruth))
  N0k <- N * (1-piY)
  N1k <- N * piY
	if(cohort == FALSE)
	{
	  N0k <- N0k * NI[1]/sum(N0k)
	  N1k <- N1k * NI[2]/sum(N1k)
	}
  popn <- as.data.frame(cbind(X, N0k, N1k))
 
  
  ## Phase I counts for the given stratification
  ##  * if strata = 'NULL' just return the population
  ##
  if(is.null(strata))
		strata <- 2:ncol(X)
	strata <- rev(sort(strata))  ## Need to reverse the ordering b/c of the way aggregate totals across the variables within 'strata'
	##
	aggCall <- paste("aggregate(popn$N0k, list(", paste("popn[,",strata,"]", sep="", collapse=", "), "), FUN=sum)", sep="")
	phaseI  <- eval(parse(text=aggCall))
	names(phaseI) <- c(names(X)[strata], "N0k")
	aggCall <- paste("aggregate(popn$N1k, list(", paste("popn[,",strata,"]", sep="", collapse=", "), "), FUN=sum)$x", sep="")
	phaseI$N1k <- eval(parse(text=aggCall))
  ##
  phaseI[,1:length(strata)] <- phaseI[,length(strata):1]  ## Reverse back again; note the colnames don't reverse so it done in the next loop
  ##
  for(i in 1:nrow(phaseI))
  	rownames(phaseI)[i] <- paste(colnames(phaseI)[length(strata):1], "=", phaseI[i, 1:length(strata)], " ", collapse = " ")
  phaseI <- round(phaseI[,-c(1:length(strata))])
	##
	cat("Expected Phase I counts:\n")
	print(phaseI)	
 
  ##
  if(!is.null(nII0) && !is.null(nII1))
  {
  	##
  	K <- prod(unlist(lapply(apply(X[,c(1,strata)], 2, unique), FUN=length)))
    if(length(nII0) != K | length(nII1) != K)
    	stop("* invalid dimension of 'nII0' or 'nII1'")
		##
		if((sum(phaseI[,1] < nII0) > 0) | (sum(phaseI[,2] < nII1) > 0))
			cat("\nWarning: At least one Phase I stratum does not have sufficiently many non-cases or cases\n\n")

		## Phase II counts
		##
		phaseII <- phaseI
		phaseII[,1] <- pmin(phaseI[,1], nII0)
		phaseII[,2] <- pmin(phaseI[,2], nII1)
		colnames(phaseII) <- c("nII0k", "nII1k")
		
		## Sampling probabilities
		##
		phaseIIprob <- phaseII/phaseI
		colnames(phaseIIprob) <- c("p0k", "p1k")		
		if(is.null(digits))
			digits <- max(-log10(phaseIIprob)) + 2
		phaseIIprob <- round(phaseIIprob, digits=digits)
		
		##
		cat("\nPhase II sample size:\n")
		print(phaseII)
		cat("\nExpected Phase II sampling probabilities:\n")
		print(phaseIIprob)
  }
  	
  ##
  invisible()
}

