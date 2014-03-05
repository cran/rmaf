\name{qn}
\alias{qn}
\title{
        Optimal and Data-Driven Moving Average Lag q
}
\description{
Determine the optimal and data-driven moving average lag \eqn{q} based on the sample size, variance and curvature of an univariate data.
}
\usage{
qn(x, data = NULL, na.action = na.fail)
}
\arguments{
  \item{x}{univariate time series or data.}
  \item{data}{an optional data frame, list or environment containing the univariate data. The default is \code{NULL}.}
  \item{na.action}{action on missing values. The default is \code{na.fail}.}
}
\details{
For an univariate series \eqn{x[t]}, the moving average filter is defined as \deqn{mhat[t] = \sum x[t]/(2q+1)} for \eqn{q + 1 \le t \le n + q}. The optimal and data-driven moving average lag \eqn{q} can be determined by using the rule-of-thumb estimator proposed in Section 3 of D. Qiu \emph{et al.} (2013). It is determined by sample size \eqn{n}, variance \eqn{\gamma(0)} and curvature \eqn{m''} of the univariate series, where \eqn{m''} is the second derivative of an unknown nonparameteric trend function. To obtain the preliminary estimators of variance \eqn{\gamma(0)} and curvature \eqn{m''}, \eqn{m} can be modeled by a cubic polynomial model. See L. Yang and R. Tscherning (1999) for more details. For the case when \eqn{q > n}, the optimal moving average lag \eqn{q} is set to be an integer part of \eqn{n^{4/5}/2}. 
}
\value{      
\item{qn}{ returns the moving average lag \eqn{q} which is used to decompose the trend, seasonality and irregularity (or residual).}
}
\references{
 D. Qiu, Q. Shao, and L. Yang (2013), Efficient inference for autoregressive coeficient in the presence of trend. \emph{Journal of Multivariate Analysis} 114, 40-53. 
 
 L. Yang, R. Tscherning (1999), Multivariate bandwidth selection for local linear regression. \emph{Journal of the Royal Statistical Society. Series B. Statistical Methodology} 61, 793-815.
}
\author{
Debin Qiu <\email{debinqiu@uga.edu}>
}
\seealso{
 \code{\link{ma.filter}}
}
\examples{
## load the global temperature data: 
## first column is time and second column is temperature.
data(globtemp)
q.n <- qn(globtemp[,2])
q.n
}
\keyword{ smooth }
\keyword{ nonparametric }
