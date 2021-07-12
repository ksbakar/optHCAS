\name{runOpt}
\alias{runOpt}
\title{Optimal pixel selection for HCAS}
\description{
This function is used to obtain co-kriging estimations based on a set of local methods. 
}
\usage{
runOpt(design.matrix, coords, response=NULL, initSampleSize = 10,
       increment.factor = 5, threshold = 0.05, optType = "exact",
       model = "univariate", weight.response = NULL, seed=1234)
}
\arguments{
\item{design.matrix}{n x p in matrix/data format.}
\item{coords}{n x 2 coordinates in matrix/data format.}
\item{response}{n x q matrix, q = 1 if univariate.}
\item{initSampleSize}{Number of initial samples.}
\item{increment.factor}{Number of increments in each iteration.}
\item{threshold}{Threshold value for convergence of the rate of change.}
\item{optType}{Can take 2 arguments: (i) "exact" for exact algorithm and (ii) "montecarlo" for montecarlo algorithm.}
\item{model}{Can take "univariate" or "multivariate".}
\item{weight.response}{Optional for multivariate model prediction, if null then use equal weights for multivariate response.}
\item{seed}{User input of seed number. Default 1234.}
}
\value{
\item{model}{Model used in the algorithm.}
\item{optType}{Algorithm type.}
\item{coords}{Pixel coordinates.}
\item{optCoords}{Optimal pixel coordinates.}
\item{optID}{Optimal ID.}
\item{comp.time}{Returns the computation time.}
}
\references{

}
%\note{
%}
\seealso{
\link{plot.pa}.
}
\examples{
\dontrun{
######################################################
## 
## simulate design matrix and coordinates
n <- 100
p <- 5
x <- c()
for(i in 1:p){
 set.seed(1234+i)
 x <- cbind(x,rnorm(n,10,1))
}
head(x)
coords <- expand.grid(x = seq(0,1,length.out=10), y = seq(0,1,length.out=10))
plot(coords,pch="*")
##
## univariate model 
beta <- c(0.2,0.5,0.4,0.3,0.6)
response <- x%*%beta + rnorm(n)
out <- runOpt(design.matrix=x, coords=coords, response=c(response),
              increment.factor = 1, threshold = 0.01, initSampleSize = 5,
              optType = "exact", model = "univariate", weight.response = NULL,
              seed=1234)
out
plot(out)
##
##
######################################################
}
}
\keyword{hcas}