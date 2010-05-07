\name{plotCcPower}
\Rdversion{1.1}
\alias{plotCcPower}
\title{
	Plot function for power, based on a case-control design
}
\description{
  The \code{plotCcPower} function plots powers of a specified
  estimate of a regression coefficient of a logistic regression model
  for different sample sizes, based on results from the function \code{ccPower} function. 

}
\usage{
plotCcPower(x, coefNum=1, yAxis=seq(from=0, to=100, by=20),
            xAxis=NULL, main=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object in a class \code{ccPower} obtained as a result of \code{\link{ccPower}} function.
}
  \item{coefNum}{
    A numeric vector number specifying the regression coefficient in \code{betaTruth} for the plot.
}
  \item{yAxis}{
    A scale marking the y-axis for the plot.
}
  \item{xAxis}{
    A scale marking the x-axis for the plot. If left as the default \code{NULL}, the x-axis scale is taken from \code{nII} in the \code{tpsResults} object
}
  \item{main}{
    Title for the plot. 
}
}
\details{
Produces a plot of statistical power (to reject a null hypothesis H0: \code{betaTruth} = 0), for the maximum likelihood estimator of a regression coefficient from a logistic regression model, based on a case-control design.
}
\author{
  Sebastien Haneuse, Takumi Saegusa
}
\seealso{
\code{\link{ccPower}}.
}
\examples{
##
data(Ohio)

##
XM   <- cbind(Int=1, Ohio[,1:3])
fitM <- glm(cbind(Death, N-Death) ~ factor(Age) + Sex + Race, data=Ohio,
            family=binomial)
betaNamesM <- c("Int", "Age1", "Age2", "Sex", "Race")

##
newBetaM      <- fitM$coef
newBetaM[2:3] <- newBetaM[2:3] / 2

\dontrun{ccResult <- ccPower(B=1000, betaTruth=newBetaM, X=XM, N=Ohio$N, r=0.5,
                    nCC=seq(from=100, to=500, by=50), betaNames=betaNamesM,
                    monitor=100)

##
par(mfrow=c(2,2))
plotCcPower(ccResult, coefNum=2, yAxis=seq(from=0, to=100, by=20),
            xAxis=seq(from=100, to=500, by=100),
            main=expression("Age effect (65-74 vs. 55-64 years), " * beta[A1]))
plotCcPower(ccResult, coefNum=3, yAxis=seq(from=0, to=100, by=20),
            xAxis=seq(from=100, to=500, by=100),
            main=expression("Age effect (75-84 vs. 55-64 years), " * beta[A2]))
plotCcPower(ccResult, coefNum=4, yAxis=seq(from=0, to=100, by=20),
            xAxis=seq(from=100, to=500, by=100),
            main=expression("Sex effect, " * beta[S]))
plotCcPower(ccResult, coefNum=5, yAxis=seq(from=0, to=100, by=20),
            xAxis=seq(from=100, to=500, by=100),
            main=expression("Race effect, " * beta[R]))
}}

