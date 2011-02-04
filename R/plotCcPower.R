plotCcPower <-
function(x,coefNum=1, yAxis=seq(from=0, to=100, by=20), xAxis=NULL, main=NULL)
{
  ##check class
  if(!class(x)=="ccPower"){
    print("Error: 'x' is not a 'ccPower' object")
    return(-1)
  }

  if(coefNum==1){
    print("Error: power for intercept is not available in a case-control study")
    invisible()
  }
  if (length(intersect(coefNum, 1:length(x$beta))) == 0) {
    print("Error: 'coefNum' is invalid")
    invisible()
    return(-1)
  }
  if(is.element(coefNum,which(x$beta==0))){
    print("Error: the specified element of 'BetaTruth' is zero")
    invisible()
    return(-1)
  }
  sub.beta <- x$beta[1:coefNum]
  num.zero <- length(sub.beta[sub.beta==0])
  coefNum <- coefNum-num.zero
  
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

