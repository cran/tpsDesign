\name{plotTpsPower}
\Rdversion{1.1}
\alias{plotTpsPower}
\title{
  Plot function for power, based on two-phase and case-control design
}
\description{
 The \code{plotTpsPower} function plots powers of a specified
 estimate of a regression coefficient of a logistic regression model for
 different sample sizes, based on results from the function \code{tpsPower} function. 
}
\usage{
plotTpsPower(x, coefNum=1, yAxis=seq(from=0, to=100, by=20),
             xAxis=NULL, main=NULL, legendXY=NULL, includeCC=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object in a class \code{tpsPower} obtained as a result of \code{\link{tpsPower}} function.
  }
  \item{coefNum}{
    A numeric vector number specifying the regression coefficient in \code{beta} for the plot.
  }
  \item{yAxis}{
    A scale marking the y-axis for the plot.
  }
  \item{xAxis}{
    A scale marking the x-axis for the plot. If left as the default \code{NULL}, the x-axis scale is taken from \code{nII} in the \code{tpsResults} object.
  }
  \item{main}{
    Title for the plot. 
  }
  \item{legendXY}{
    Optional vector indicating the co-ordinates for the top-left hand corner of the legend box.
  }
  \item{includeCC}{
    Logical indicator of whether or not to include power for the case-control (CC) estimator.
}
}
\details{
Produces a plot of statistical power (to reject a null hypothesis H0: beta = 0), for estimators of a regression coefficient from a logistic regression model, based on a case-control design.
}
\author{
  Sebastien Haneuse, Takumi Saegusa
}
\seealso{
  \code{\link{tpsPower}}.
}
\examples{
##
data(Ohio)

##
XM   <- cbind(Int=1, Ohio[,1:3])
fitM <- glm(cbind(Death, N-Death) ~ factor(Age) + Sex + Race, data=Ohio,
            family=binomial)
betaNamesM <- c("Int", "Age1", "Age2", "Sex", "Race")

## Power for the TPS design where phase I stratification is based on Age
newBetaM      <- fitM$coef
newBetaM[2:3] <- newBetaM[2:3] / 2

\dontrun{powerAgeTPS <- tpsPower(B=1000, betaTruth=newBetaM, X=XM, N=Ohio$N,
                        strata=2, nII=seq(from=100, to=500, by=50),
                        betaNames=betaNamesM, monitor=100)
plotTpsPower(powerAgeTPS, coefNum=3, yAxis=seq(from=0, to=100, by=20),
             xAxis=seq(from=100, to=500, by=100),
             main=expression("Age effect (75-84 vs. 55-64 years), "
             * beta[A2]), legendXY=c(400, 70))
}
## Power for the TPS design where phase I stratification is based on Race.

\dontrun{powerRaceTPS <- tpsPower(B=1000, betaTruth=fitM$coef, X=XM, N=Ohio$N,
                         strata=4, nII=seq(from=100, to=1000, by=100),
                         betaNames=betaNamesM, monitor=100)
par(mfrow=c(2,2))
plotTpsPower(powerRaceTPS, coefNum=2, yAxis=seq(from=0, to=100, by=20),
             xAxis=seq(from=100, to=1000, by=100),
             main=expression("Age effect (65-74 vs. 55-64 years), "
             * beta[A1]), legendXY=c(800, 60))
plotTpsPower(powerRaceTPS, coefNum=3, yAxis=seq(from=0, to=100, by=20),
             xAxis=seq(from=100, to=1000, by=100),
             main=expression("Age effect (75-84 vs. 55-64 years), "
             * beta[A2]), legendXY=c(800, 60))
plotTpsPower(powerRaceTPS, coefNum=4, yAxis=seq(from=0, to=100, by=20),
             xAxis=seq(from=100, to=1000, by=100),
             main=expression("Sex effect, " * beta[S]),
             legendXY=c(800, 60))
plotTpsPower(powerRaceTPS, coefNum=5, yAxis=seq(from=0, to=100, by=20),
             xAxis=seq(from=100, to=1000, by=100),
             main=expression("Race effect, " * beta[R]),
             legendXY=c(800, 60))
}
}


