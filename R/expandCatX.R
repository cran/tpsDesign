expandCatX <-
function(X)
{
	## Assumes a check has been performed to make sure the columns adhere to the {0,1,2,...} coding convention
	##
	value <- matrix(X[,1], ncol=1, dimnames=list(1:12, "Int"))
	##
	if(ncol(X) > 1)
	{
		n.lev <- unlist(lapply(apply(X, 2, unique), FUN=length))
		for(i in 2:ncol(X))
		{
			for(j in 1:(n.lev[i]-1))
			{
				value <- cbind(value, as.numeric(X[,i] == j))
				if(n.lev[i] == 2) colnames(value)[ncol(value)] <- colnames(X)[i]
				if(n.lev[i] > 2) colnames(value)[ncol(value)] <- paste(colnames(X)[i], ".", j, sep="")
			}
		}
	}
	##
	return(value)
}

