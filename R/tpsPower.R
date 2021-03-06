tpsPower <-
function(B=1000,
										 betaTruth,
										 X,
										 N,
										 strata,
										 nII,
										 alpha=0.05,
                     digits=1,
                     betaNames=NULL,
                     monitor=NULL,
                     cohort=TRUE,
                     NI=NULL)
{
	## Checks
	##
	problem <- coreChecks(betaTruth=betaTruth, X=X, N=N, betaNames=betaNames)
	if(problem != "")
		stop(problem)
 	##
	problem <- tpsChecks(X=X, strata=strata, nII=nII, cohort=cohort, NI=NI)
	if(problem != "")
		stop(problem)
 	##
  if(is.null(colnames(X)))
  	colnames(X) <- c("Int", paste("V", 1:(ncol(X) - 1), sep = ""))
	if(is.null(monitor))
		monitor <- B + 1

	## Expanded design matrix, outcome probabilities and formulae for glm() and tps() calls
  ##
  Xexp <- as.data.frame(expandCatX(X))
  if(!is.null(betaNames))
  	names(Xexp) <- betaNames
  piY <- expit(as.vector(as.matrix(Xexp) %*% betaTruth))
	formCD  <- as.formula(paste("cbind(N1, N0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))
	formTPS <- as.formula(paste("cbind(n1, n0) ~", paste(colnames(Xexp)[-1], collapse=" + ", sep="")))

	## Phase I stratification variable
	##
	Xexp$S <- stratify(X, strata)
	K      <- length(unique(Xexp$S))

	## Run simulation
	##
	lenII    <- length(nII)
  nDesigns <- 1 + (4*lenII)                                      ## CD plus (CC, WL, PL, ML) for each value of NII
  p        <- length(betaTruth)
	waldTest <- array(NA, dim=c(B, nDesigns, p))
  cat(paste(1+(2*lenII), "designs will be simulated\n"))
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
		waldTest[b,1,] <- (fitCD[,4] < alpha)
		
		
		##
   	op <- options()
		##
		for(i in 1:lenII)
		{
			## Case-control study design
	    ##
  		n0 <- round(nII[i] * 0.5)
  		n1 <- round(nII[i] * 0.5)
			##
			Xexp$n0 <- rmvhyper(Xexp$N0, n0)
			Xexp$n1 <- rmvhyper(Xexp$N1, n1)
			##
			index <- 1 + i + (0*lenII)
			fitCC <- try(glm(formTPS, data=Xexp, family=binomial()), silent=TRUE)
			if(class(fitCC)[1] == "glm")
			{
				fitCC <- summary(fitCC)$coef
				if(nrow(fitCC) == p)
					waldTest[b,index,] <- c(FALSE, (fitCC[-1,4] < alpha))
			}

			## Two-phase study design
			##
	    ## Phase I counts
  	  ##
    	if(cohort == FALSE)
    	{
				Xexp$N0 <- rmvhyper(Xexp$N0, NI[1])
				Xexp$N1 <- rmvhyper(Xexp$N1, NI[2])
    	}
   	 	nn0 <- tapply(Xexp$N0, Xexp$S, FUN=sum)
   		nn1 <- tapply(Xexp$N1, Xexp$S, FUN=sum)
			##
  		n0 <- round(nII[i] * 0.5 / K)
  		n1 <- round(nII[i] * 0.5 / K)
	    ##
	    Xexp$n0 <- Xexp$N0
	    Xexp$n1 <- Xexp$N1
	    for(k in 1:K)
	    {
  	    Xexp$n0[Xexp$S == k] <- rmvhyper(Xexp$N0[Xexp$S == k], n0)
  	    Xexp$n1[Xexp$S == k] <- rmvhyper(Xexp$N1[Xexp$S == k], n1)
    	}
			##
			index <- 1 + i + (1*lenII)
    	options(warn=-1)
    	fitWL <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="WL", cohort=cohort), silent=TRUE)
    	options(op)
			if(class(fitWL)[1] == "tps")
				waldTest[b,index,] <- abs(fitWL$coef/sqrt(diag(fitWL$cove))) > abs(qnorm(alpha/2))
    	##
			index <- 1 + i + (2*lenII)
    	fitPL <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="PL", cohort=cohort), silent=TRUE)
			if(class(fitPL)[1] == "tps")
				waldTest[b,index,] <- abs(fitPL$coef/sqrt(diag(fitPL$covm))) > abs(qnorm(alpha/2))
    	##
			index <- 1 + i + (3*lenII)
    	fitML <- try(tps(formTPS, Xexp, nn0=nn0, nn1=nn1, Xexp$S, method="ML", cohort=cohort), silent=TRUE)
			if(class(fitML)[1] == "tps")
				waldTest[b,index,] <- abs(fitML$coef/sqrt(diag(fitML$covm))) > abs(qnorm(alpha/2))
		}
	}
	
	## Error checks for output and evaluate power
	##
	if(lenII == 1)
		colNames <- c("CD  ", "CC  ", "TPS-WL  ", "TPS-PL  ", "TPS-ML  ")
	if(lenII > 1)
	{
		colNames <- "CD"
		colNames <- c(colNames, paste("CC     nII =", nII[1], " "))
		for(i in 2:lenII) colNames <- c(colNames, paste("       nII =", nII[i], " "))
		colNames <- c(colNames, paste("TPS-WL nII =", nII[1], " "))
		for(i in 2:lenII) colNames <- c(colNames, paste("       nII =", nII[i], " "))
		colNames <- c(colNames, paste("TPS-PL nII =", nII[1], " "))
		for(i in 2:lenII) colNames <- c(colNames, paste("       nII =", nII[i], " "))
		colNames <- c(colNames, paste("TPS-ML nII =", nII[1], " "))
		for(i in 2:lenII) colNames <- c(colNames, paste("       nII =", nII[i], " "))
	}
	##
	keep    <- apply(is.na(waldTest), c(1,2), sum) == 0
	results <- evalPower(waldTest, keep, resNames=list(colNames, betaNames))
	results[results[,1] == 0, 1] <- NA
	failed <- B - matrix(apply(keep, 2, sum), ncol=1, dimnames=list(rownames(results), ""))
	
	## Return object of class 'tpsPower'
  ##
  value           <- NULL
  value$B         <- B
  value$betaTruth <- betaTruth
  value$X         <- X
  value$N         <- N
  value$strata    <- strata
  value$nII       <- nII
  value$alpha     <- alpha
	value$digits 	  <- digits
	value$cohort    <- cohort
	value$NI        <- NI
  ##
  value$failed    <- failed
  value$betaPower <- results
  ##
  class(value) <- "tpsPower"
  return(value)
}

