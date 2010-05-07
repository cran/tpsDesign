plotTpsPower <-
function(x,coefNum=1, yAxis=seq(from=0, to=100, by=20), xAxis=NULL, main=NULL, legendXY=NULL, includeCC=TRUE)
{
  ##check class
  if(!class(x)=="tpsPower"){
    print("Error: 'x' is not a 'tpsPower' object")
    return(0)
  }
  
  ##check coefNum
  if(!is.element(coefNum,1:length(x$beta))){
    print("Error: 'coefNum' is invalid")
    return(-1)
  }
   
  ##
  if(is.null(xAxis)) xAxis <- x$nII
  if(is.null(title)) title <- paste("Power for", colnames(x$power)[coefNum])
  
  ##
  nLvls   <- length(x$nII)
  ##
  if(nLvls>1){
    powerCC  <- x$power[(2 + c(0:(nLvls-1))*4), coefNum]
    powerWL  <- x$power[(3 + c(0:(nLvls-1))*4), coefNum]
    powerPL  <- x$power[(4 + c(0:(nLvls-1))*4), coefNum]
    powerML  <- x$power[(5 + c(0:(nLvls-1))*4), coefNum]
  }else{
    powerCC <- x$power[2, coefNum]
    powerWL <- x$power[3, coefNum]
    powerPL <- x$power[4, coefNum]
    powerML <- x$power[5, coefNum]
  }   
  
  ##
  plot(range(xAxis), range(yAxis), xlab="Phase II sample size, n", ylab="Power", main=main, type="n", axes=FALSE)
  axis(1, at=xAxis)
  axis(2, at=yAxis)
  points(x$nII, powerWL, pch=2)
  points(x$nII, powerPL, pch=3)
  points(x$nII, powerML, pch=4)
  lines(x$nII, powerWL, lty=1, lwd=2)
  lines(x$nII, powerPL, lty=2, lwd=2)
  lines(x$nII, powerML, lty=3, lwd=2)
  ##
  if(includeCC == TRUE) lines(x$nII, powerCC, lty=4, lwd=2)
  ##
  if(is.null(legendXY)){
    if(includeCC==FALSE){
      legendXY <- c(max(xAxis) - ((max(xAxis)-min(xAxis)) * 0.2), min(yAxis) + ((max(yAxis)-min(yAxis))*0.18))
    }else{
      legendXY <- c(max(xAxis) - ((max(xAxis)-min(xAxis)) * 0.2), min(yAxis) + ((max(yAxis)-min(yAxis))*0.2))
    }
  }
  if(includeCC == FALSE) legend(legendXY[1], legendXY[2], c("WL", "PL", "ML"), pch=c(2:4), lwd=c(2,2,2), lty=c(1:3))
  if(includeCC == TRUE) legend(legendXY[1], legendXY[2], c("WL", "PL", "ML", "CC"), pch=c(2:4,NA), lwd=c(2,2,2,2), lty=c(1:4))
  ##
  
  invisible()
}

