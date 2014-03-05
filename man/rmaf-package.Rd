\name{rmaf-package}
\alias{rmaf-package}
\alias{rmaf}
\docType{package}
\title{
Refined Moving Average Filter Package
}
\description{
A refined moving average filter using the optimal and data-driven moving average lag \eqn{q} to decompose the trend, seasonality and irregularity for an univariate time series or data.
}
\details{
\tabular{ll}{
Package: \tab rmaf\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-03-05\cr
License: \tab GPL (>= 2)\cr
}
This package contains a function to determine the optimal and data-driven moving average lag \eqn{q}, and function to decompose the trend, seasonality and irregularity for an univariate time series.  A dataset of the first differences of annual global surface air temperatures in Celsius from 1880 through 1985 is also included in the package for illustrating the trend decomposition. 

For a complete list of functions and dataset, use \code{library(help = rmaf)}.
}
\author{
Debin Qiu

Maintainer: Debin Qiu <\email{debinqiu@uga.edu}>

}
\references{
D. Qiu, Q. Shao, and L. Yang (2013), Efficient inference for autoregressive coeficient in the presence of trend. \emph{Journal of Multivariate Analysis} 114, 40-53. 

P.J. Brockwell, R.A. Davis, Time Series: Theory and Methods, second ed., Springer, New York, 1991.
}

\seealso{
\code{\link{ma.filter}}, \code{\link{qn}}}
\keyword{ smooth }
\examples{
## The first difference of annual global surface air temperatures from 1880 to 1985 with only trend
data(globtemp)
x <- globtemp[,2]
q.n <- qn(x)
fit <- ma.filter(x)
}
