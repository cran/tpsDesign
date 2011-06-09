beta0 <-
function(betaX,
									X,
									N,
									rhoY)
{
	##
  value <- optimize(beta0eval,
  									interval=logit(rhoY) + c(-2,2),
                    betaX=betaX,
                    X=X,
                    N=N,
                    rhoY=rhoY,
                    tol=.Machine$double.eps^0.5)$minimum
  ##
  return(value)
}

