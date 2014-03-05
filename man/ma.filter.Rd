\name{ma.filter}
\alias{ma.filter}
\title{
Refined Moving Average Filter to Decompose Trend, Seasonality and Irregularity
}
\description{
Use refined moving average filter to decompose the trend, seasonality and irregularity for an univariate time series. 
}
\usage{
ma.filter(xt, data = NULL, seasonal = FALSE, d = NULL, na.action = na.fail, plot = FALSE)
}
\arguments{
  \item{xt}{
            univariate time series to be decomposed.
}
  \item{data}{
            an optional data frame, list or environment containing the univariate time series. The default is \code{NULL}.
}
\item{seasonal}{
            seasonality. The default is \code{FALSE} which is to decompose the trend only. When \code{seasonal} is \code{TRUE}, the seasonal period \code{d} should be specified and an additive seasonality is only considered here. }
  \item{d}{
            seasonal period. It should be at least 2 and divisor of sample size. The default is \code{NULL}.
}
  \item{na.action}{
            action on missing values. The detault is \code{na.fail}.
}
  \item{plot}{
            plots of data, trend, seasonality (if seasonal \eqn{d} is specified) and irregularity. The default is \code{FALSE}. 
}
}
\details{
For an univariate time seties \eqn{x[t]}, the additive seasonal model is assumed to be \deqn{x[t] = m[t] + S[t] + R[t],} where \eqn{m[t], S[t], R[t]} are trend, seasonality and irregular components, respectively. The trend and seasonality of an univariate time series \eqn{x[t]} are decomposed by using the refined moving average filter. Once we obtain the optimal moving average lag \eqn{q}, the trend can be estimated by using the refined moving average \deqn{mhat[t] = \sum x[t]/(2q+1),} for \eqn{q + 1 \le t \le n - q}. If \eqn{q + 1 > n - q}, we take \eqn{q = min(n - q, q)}.
 If there is no seasonal component, the irregularity or residual can be obtained by \eqn{Rhat[t] = x[t] - mhat[t]}. Otherwise, the seasonality \eqn{Shat[t]} can be decomposed by using the refined moving average on \eqn{x[t] - mhat[t]} and then irregularity or residual is \eqn{Rhat[t] = x[t] - mhat[t] - Shat[t]}. For \eqn{t < q + 1} and \eqn{t > n - q}, the corresponding estimators are based on equation (7) in D. Qiu \emph{et al.} (2013). More details about decomposing trend and seasonality can be seen in Section 1.5 of P.J. Brockwell \emph{et al.} (1991).

For the multiplicative seasonal model \deqn{x[t] = m[t] * S[t] * R[t],} it can be transformed to an additive seasonal model by taking the logarithm on both sides if \eqn{x[t]>0}, i.e., \deqn{log(x[t]) = log(m[t]) + log(S[t]) + log(R[t]),} and then use the refined moving average filter for the components decomposition as in the additive seasonal model. 
}
\value{
A list of time series decomposition components including:
  \item{data }{original univariate time series.}
  \item{trend }{fitted trend.}
  \item{season}{fitted seasonality.}
  \item{remainder}{fitted irregularity or residual.}
  \item{SSE}{sum of squared errors.}
}
\references{
 D. Qiu, Q. Shao, and L. Yang (2013), Efficient inference for autoregressive coeficient in the presence of trend. \emph{Journal of Multivariate Analysis} 114, 40-53. 
 
 P.J. Brockwell, R.A. Davis, Time Series: Theory and Methods, second ed., Springer, New York, 1991.
}
\author{
Debin Qiu <\email{debinqiu@uga.edu}>
}

\seealso{
\code{\link{qn}}
}
\examples{
## decompose the trend for the first difference of annual global air temperature from 1880-1985
data(globtemp)
decomp1 <- ma.filter(globtemp[,2], plot = TRUE)

## decompose the trend and seasonality for CO2 data with monthly and additive seasonality
data(co2)
decomp2 <- ma.filter(co2, d = 12, seasonal = TRUE, plot = TRUE)

## decompose the trend and seasonality for monthly airline passenger numbers from 1949-1960
data(AirPassengers)
decomp3 <- ma.filter(AirPassengers, d = 12, seasonal = TRUE, plot = TRUE)

## simulation data: oracally efficient estimation for AR(p) coefficients
d <- 12
n <- d*100
x <- (1:n)/n
y <- 1 + 2*x + 0.3*x^2 + sin(pi*x/6) + arima.sim(n = n,list(ar = 0.2), sd = 1)
fit <- ma.filter(y, d = 12, seasonal = TRUE)
ar(fit$remainder, aic = FALSE, order.max = 1)$ar
}
\keyword{ smooth }
\keyword{ nonparametric }