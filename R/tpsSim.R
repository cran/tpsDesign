tpsSim <-
function(B=1000,
									 betaTruth,
									 X,
									 N,
									 strata,
									 nII0=NULL,
									 nII1=NULL,
									 nII=NULL,
                   nCC=NULL,
                   alpha=0.05,
                   threshold=c(-Inf,Inf),
                   digits=1,
                   betaNames=NULL,
                   referent=2,
                   monitor=NULL,
                   cohort=TRUE,
                   NI=NULL,
                   returnRaw=FALSE)
{
	## Checks
	##
	problem <- coreChecks(betaTruth=betaTruth, X=X, N=N, betaNames=betaNames)
	if(problem != "")
		stop(problem)
	##
	problem <- tpsChecks(X=X, strata=strata, nII=nII, cohort=cohort, NI=NI, nII0=nII0, nII1=nII1, nCC=nCC, threshold=threshold)
	if(problem != "")
		stop(problem)
 	##
  if(is.null(colnames(X)))
  	colnames(X) <- c("Int", paste("V", 1:(ncol(X) - 1), sep = ""))
	if(is.null(monitor))
		monitor <- B + 1
	if(is.null(nCC))
	{
		if(!is.null(nII))
			nCC <- nII
		else
			nCC <- c(sum(nII0), sum(nII1))
	}

	## Expanded design matrix, outcome probabilities and formulae for glm() and tps() calls
  ##
  Xexp <- as.data.frame(expandCatX(X))
  if(!is.null(betaNames))
  	names(Xexp) <- betaNames
  piY <- expit(as.vector(as.matrix(Xexp) %*% betaTruth))
	formCD  <- as.formula(paste("cbind(N1, N0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))
	formTPS <- as.formula(paste("cbind(n1, n0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))

	## Phase I stratification(s)
	##
	if(max(strata) > 0)
	{
		strataMat <- as.matrix(stratify(X, strata), ncol=1)
	}
	if(max(strata) == 0)
	{
		colIndex  <- 2:ncol(X)
		strataLab <- NULL
		strataMat <- NULL
		for(i in 1:(length(colIndex)-1))
		{
			temp <- combn(colIndex, i)
			for(j in 1:ncol(temp))
			{
				vec10     <- 10^((nrow(temp)-1):0)
				strataLab <- c(strataLab, paste(sum(temp[,j]*vec10), paste(rep(" ", length(colIndex)-1-i), collapse=""), sep=""))
				strataMat <- cbind(strataMat, stratify(X, strata=temp[,j]))
			}
		}
	}
	##
	nStrat <- ncol(strataMat)
	
	## Run simulation
	##
  nDesigns <- 1 + 1 + (3*nStrat)                      ## CD plus CC plus (WL, PL, ML) for each two-phase design
  p        <- length(betaTruth)
	betaHat  <- array(NA, dim=c(B, nDesigns, p))
	seHat    <- array(NA, dim=c(B, nDesigns, p))
	waldTest <- array(NA, dim=c(B, nDesigns, p))
  cat(paste(2+nStrat, "designs will be simulated\n"))
	##
	for(b in 1:B)
	{
		##
    if((b %% monitor) == 0)
      cat("Repetition", b, "of", B, "complete\n")
    
		## Complete data study design
		##
		Xexp$N1 <- rbinom(nrow(Xexp), N, piY)
		Xexp$N0 <- N - Xexp$N1
		fitCD   <- summary(glm(formCD, data=Xexp, family=binomial()))$coef
		betaHat[b,1,]  <- fitCD[,1]
		seHat[b,1,]    <- fitCD[,2]
		waldTest[b,1,] <- (fitCD[,4] < alpha)
		
		## Case-control study design
		##
		Xexp$n0 <- rmvhyper(Xexp$N0, nCC[1])
		Xexp$n1 <- rmvhyper(Xexp$N1, nCC[2])
		##
		fitCC <- try(glm(formTPS, data=Xexp, family=binomial()), silent=TRUE)
		if(class(fitCC)[1] == "glm")
		{
			fitCC <- summary(fitCC)$coef
			if(nrow(fitCC) == p)
			{
				betaHat[b,2,-1] <- fitCC[-1,1]
				seHat[b,2,-1]   <- fitCC[-1,2]
				waldTest[b,2,-1] <- (fitCC[-1,4] < alpha)
			}
		}

		## Two-phase study design(s)
		##
		##
   	op <- options()
		##
		for(s in 1:nStrat)
		{
			##
			Xexp$S <- strataMat[,s]
			K      <- length(unique(Xexp$S))
			
    	## Phase I counts
    	##
    	if(cohort == FALSE)
    	{
				Xexp$N0 <- rmvhyper(Xexp$N0, NI[1])
				Xexp$N1 <- rmvhyper(Xexp$N1, NI[2])
    	}
  	  nn0 <- tapply(Xexp$N0, Xexp$S, FUN=sum)
	    nn1 <- tapply(Xexp$N1, Xexp$S, FUN=sum)

			## Phase II data
	   	##
	    Xexp$n0 <- Xexp$N0
	    Xexp$n1 <- Xexp$N1
	    for(k in 1:K)
	    {
  	    ## Controls
  	    if(is.null(nII0))
  	    	Xexp$n0[Xexp$S == k] <- rmvhyper(Xexp$N0[Xexp$S == k], round(nII[1]/K))
  	    if(!is.null(nII0))
  	    	Xexp$n0[Xexp$S == k] <- rmvhyper(Xexp$N0[Xexp$S == k], nII0[k])
  	  	## Cases
  	    if(is.null(nII1))
  	    	Xexp$n1[Xexp$S == k] <- rmvhyper(Xexp$N1[Xexp$S == k], round(nII[2]/K))
  	    if(!is.null(nII1))
  	    Xexp$n1[Xexp$S == k] <- rmvhyper(Xexp$N1[Xexp$S == k], nII1[k])
    	}
    	
    	## Estimation and inference
			##
			index <- 2 + s + (0*nStrat)
    	options(warn=-1)
    	fitWL <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="WL", cohort=cohort), silent=TRUE)
    	options(op)
			if(class(fitWL)[1] == "tps")
			{
				betaHat[b,index,]  <- fitWL$coef
				seHat[b,index,]    <- sqrt(diag(fitWL$cove))
				waldTest[b,index,] <- abs(fitWL$coef/sqrt(diag(fitWL$cove))) > abs(qnorm(alpha/2))
			}
    	##
			index <- 2 + s + (1*nStrat)
    	fitPL <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="PL", cohort=cohort), silent=TRUE)
			if(class(fitPL)[1] == "tps")
			{
				betaHat[b,index,]  <- fitPL$coef
				seHat[b,index,]    <- sqrt(diag(fitPL$covm))
				waldTest[b,index,] <- abs(fitPL$coef/sqrt(diag(fitPL$covm))) > abs(qnorm(alpha/2))
			}
    	##
			index <- 2 + s + (2*nStrat)
    	fitML <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="ML", cohort=cohort), silent=TRUE)
			if(class(fitML)[1] == "tps")
			{
				betaHat[b,index,]  <- fitML$coef
				seHat[b,index,]    <- sqrt(diag(fitML$covm))
				waldTest[b,index,] <- abs(fitML$coef/sqrt(diag(fitML$covm))) > abs(qnorm(alpha/2))
			}
		}
	}
	
	## Error checks for output and evaluate operating characteristics
	##
	if(nStrat == 1)
		colNames <- c("CD ", "CC ", "TPS-WL ", "TPS-PL ", "TPS-ML ")
	if(nStrat > 1)
	{
		colNames <- c("CD", "CC")
		colNames <- c(colNames, paste("TPS-WL", strataLab[1], ""))
		for(s in 2:nStrat) colNames <- c(colNames, paste("      ", strataLab[s], ""))
		colNames <- c(colNames, paste("TPS-PL", strataLab[1], ""))
		for(s in 2:nStrat) colNames <- c(colNames, paste("      ", strataLab[s], ""))
		colNames <- c(colNames, paste("TPS-ML", strataLab[1], ""))
		for(s in 2:nStrat) colNames <- c(colNames, paste("      ", strataLab[s], ""))
	}
	keep     <- keepOC(betaTruth, betaHat, seHat, threshold)
	results  <- evalOC(betaTruth, betaHat, seHat, waldTest, keep, alpha=alpha, referent=referent, resNames=list(colNames, betaNames))
	failed   <- B - matrix(apply(keep, 2, sum), ncol=1, dimnames=list(rownames(results$betaPower), ""))
	
	## Return object of class 'tpsPower'
  ##
  value           <- NULL
  value$B         <- B
  value$betaTruth <- betaTruth
  value$X         <- X
  value$N         <- N
  value$strata    <- strata
  value$nII0      <- nII0
  value$nII1      <- nII1
  value$nII       <- nII
  value$nCC       <- nCC
  value$alpha     <- alpha
  value$threshold <- threshold
	value$digits 	  <- digits
	value$cohort    <- cohort
	value$NI        <- NI
  ##
  value$failed    <- failed
  value$results   <- results
  ##
  if(returnRaw == TRUE)
  {
  	value$betaHat <- betaHat
		value$seHat   <- seHat
  }
  class(value) <- "tpsSim"
  return(value)
}

