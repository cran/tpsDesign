plotCcPower <-
function(x,coefNum=1, yAxis=seq(from=0, to=100, by=20), xAxis=NULL, main=NULL)
{
  ##check class
  if(!class(x)=="ccPower"){
    print("Error: 'x' is not a 'ccPower' object")
    return(0)
  }

  ##
  if(is.null(xAxis)) xAxis <- x$nII
  if(is.null(title)) title <- paste("Power for", colnames(x$power)[coefNum])
  
  ##
  powerCC <- x$power[-1, coefNum]
  
  ##
  plot(range(xAxis), range(yAxis), xlab="Case-control sample size, n", ylab="Power", main=main, type="n", axes=FALSE)
  axis(1, at=xAxis)
  axis(2, at=yAxis)
  points(x$nII, powerCC, pch=2)
  lines(x$nII, powerCC, lty=1, lwd=2)
  ##
  invisible()
}

