\name{plot.pa}
\alias{plot.pa}
\title{Summary statistics of the parameters.}
\description{
This function is used to obtain MCMC summary statistics.
}
\usage{
\method{plot}{pa}(x, treatment=FALSE, differences=FALSE, zval=FALSE, 
        variance=FALSE, covariance=FALSE, length.out=10, ...)
}
\arguments{
\item{x}{Object of class inheriting from "pa".}
\item{treatment}{Logical, if TRUE then plots the treatment surface.}
\item{differences}{Logical, if TRUE then plots the surface of the treatment differences.}
\item{zval}{Logical, if TRUE then plots the surface of the Z-value for the treatment differences.}
\item{variance}{Logical, if TRUE then plots the surface of the treatment variances.}
\item{covariance}{Logical, if TRUE then plots the surface of the treatment covariances.}
\item{length.out}{Numeric statement to adjust the labels of the plot. Default is 10.}
\item{...}{Other arguments.}
}
\value{
}
\seealso{
\code{\link{run}}.
}
\examples{
\dontrun{
##

plot(out, treatment=TRUE) # where out is the output from pa class
plot(out, difference=TRUE) # 
plot(out, zval=TRUE) # 
plot(out, variance=TRUE) # 
plot(out, covariance=TRUE) # 

##
}
}
\keyword{pa}
