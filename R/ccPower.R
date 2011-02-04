ccPower <-
function (B = 1000, betaTruth, X, N, nCC, r,
                     alpha = 0.05,threshold = c(-Inf, Inf),
                     digits = NULL, betaNames = NULL,
                     monitor = NULL, NI=NULL){
  nII <- nCC
  beta <- betaTruth
  if(min(N, nII) < 0) {
    print("Error: sample size is not positive")
    return(-1)
  }
  if (min(X) < 0) {
    print("Error: invalid design matrix 'X'")
    return(-1)
  }
  if(length(N) != nrow(X)) {
    print("Error: invalid dimensions of 'N' and 'X'")
    return(-1)
  }
  if(length(r) != 1) {
    print("Error: 'r' is not a single number")
    return(-1)
  }
  nlevs <- 1
  for(i in 1:(ncol(X) - 1)) {
    nlevs <- c(nlevs, length(unique(X[, (i + 1)])) - 1)
  }
  if(length(beta) != sum(nlevs)) {
    print("Error: invalid dimension of 'beta'")
    return(-1)
  }
  if(!is.null(threshold[1])) {
    if(length(threshold) != 2) {
      print("Error: 'threshold' is not a pair of numbers")
    }
  }
  if(!is.null(betaNames[1])) {
    if(length(beta) != length(betaNames)) {
      print("Error: 'beta' and 'betaNames' are not of the same length")
      return(-1)
    }
    if(!is.character(betaNames)) {
      print("Error: elements of 'betaNames' are not character")
      return(-1)
    }
  }
  ## case-control option
  if(!is.null(NI[1])){
    if(length(NI)!=2){
      print("Error: 'NI' is a pair of Phase I sample sizes for controls and cases")
      return(-1)
    }
    if(min(NI)<0){
      print("Error: sample size is not positive")
      return(-1)
    }
  } 
     
  if(!is.null(colnames(X))) {
    var.name <- colnames(X)[-1]
    if(is.null(betaNames)) {
      if(is.null(names(beta))) {
        betaNames <- "Intercept"
        for(i in 2:ncol(X)) {
          betaNames <- c(betaNames, paste("(", var.name[i - 1], ".", 1:nlevs[i], ")", sep = ""))
        }
        names(beta) <- betaNames
      }else{
        betaNames <- names(beta)
      }
    }
  }else{
    colnames(X) <- c("Intercept", paste("V", 1:(ncol(X) - 1), sep = ""))
    var.name <- colnames(X)[-1]
    if (is.null(betaNames)) {
      if (is.null(names(beta))) {
        betaNames <- "Intercept"
        var.name2 <- paste("(V", 1:(ncol(X) - 1), ")", sep = "")
        for (i in 2:ncol(X)) {
          betaNames <- c(betaNames, paste(var.name[i - 1], ".", 1:nlevs[i], sep = ""))
        }
        names(beta) <- betaNames
      }else{
        betaNames <- names(beta)
      }
    }
  }
  if(is.null(monitor)) {
    monitor <- B + 1
  }else{
    monitor <- as.integer(monitor)
    if (monitor < 0) {
      print("Error: 'monitor' is not a positive number")
      return(-1)
    }
  }
  num.nii <- length(nII)
  result <- NULL
  result1 <- NULL
  result2 <- NULL
  beta.0.ind <- which(beta == 0)
  n0 <- round(nII[1] * r/(r + 1))
  n1 <- nII[1] - n0
  print(paste("There are", num.nii, "case-control designs to simulate"))
  result1 <- ccSimOne(B = B, betaTruth = beta, X = X, N = N, 
                      n0 = n0, n1 = n1, alpha = alpha, threshold = threshold, 
                      monitor = monitor,NI=NI)
  print("Design 1 complete")
  if (num.nii > 1) {
    n0 <- round(nII[2] * r/(r + 1))
    n1 <- nII[2] - n0
    result2 <- ccSimOne(B = B, betaTruth = beta, X = X, N = N, 
                        n0 = n0, n1 = n1, alpha = alpha, threshold = threshold, 
                        monitor = monitor,NI=NI)
    print("Design 2 complete")
    result <- list(result1, result2)
    if (num.nii > 2) {
      for (i in 3:num.nii) {
        n0 <- round(nII[i] * r/(r + 1))
        n1 <- nII[i] - n0
        result[[i]] <- ccSimOne(B = B, betaTruth = beta, 
                                X = X, N = N, n0 = n0, n1 = n1, alpha = alpha, 
                                threshold = threshold, monitor = monitor,NI=NI)
        print(paste("Design", i, "complete."))
      }
    }
  }else{
    result <- result1
  }
  colname <- betaNames
  rowname <- c("CD", rep(c(" CC"), times = num.nii))
  max.numnii <- max(nchar(nII))
  spc <- NULL
  for(i in 1:max.numnii) {
    spc <- paste(spc, " ", sep = "")
  }
  for(i in 1:num.nii) {
    for (j in 1:(max.numnii - 1)) {
      if (nchar(nII[i]) == j) {
        for (k in 1:(max.numnii - j)) {
          nII[i] <- paste(nII[i], " ", sep = "")
        }
        break
      }
    }
  }
  for(i in 1:num.nii) {
    rowname[(i + 1)] <- paste(nII[i], rowname[(i + 1)], sep = "")
  }
  matrix.power <- NULL
  if(num.nii != 1) {
    matrix.power <- result[[1]]$power[1:2, ]
    for (i in 2:num.nii) {
      matrix.power <- rbind(matrix.power, result[[i]]$power[2, 
                                                            ])
    }
  }else{
    matrix.power <- result$power[1:2, ]
  }
  rownames(matrix.power) <- rowname
  if(length(beta.0.ind) == 0) {
    colnames(matrix.power) <- colname
  }else{
    colnames(matrix.power) <- colname[-beta.0.ind]
  }
  matrix.na <- NULL
  if (num.nii != 1) {
    matrix.na <- result[[1]]$na[, 1:2]
    for(i in 2:num.nii) {
      matrix.na <- cbind(matrix.na, result[[i]]$na[, 2])
    }
  }else{
    matrix.na <- result$na[, 1:2]
  }
  matrix.na <- t(matrix.na)
  colnames(matrix.na) <- c("Point Estimate", "Above Threshold", 
                           "Standard Deviation")
  rownames(matrix.na) <- rowname
  output <- NULL
  output$B <- B
  output$beta <- beta
  output$X <- X
  output$N <- N
  output$nII <- as.numeric(nII)
  output$r <- r
  output$alpha <- alpha
  output$threshold <- threshold
  if(is.null(digits)) {
    output$digits <- 0
  }else{
    output$digits <- digits
  }
  if(!is.null(NI[1])){
    output$NI <- NI
  }
  output$power <- matrix.power
  output$na <- matrix.na
  class(output) <- "ccPower"
  return(output)
}

