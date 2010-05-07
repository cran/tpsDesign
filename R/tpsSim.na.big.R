tpsSim.na.big <-
function(x,betaTruth,threshold){
  big <- 0
  if(sum(is.na(x))>0)
  return(0)
  big <- ifelse(sum(x-betaTruth<threshold[1]|x-betaTruth>threshold[2])>0,1,0)
  return(big)
}

