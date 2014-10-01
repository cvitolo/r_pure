\name{gisInfo}
\alias{gisInfo}
\docType{data}
\title{
gisInfo
}
\description{
A data.frame that describes gauges on a catchment. It can be used to test the function ArealInterpolation with the methods "IDW" and "OK".
}
\usage{data(gisInfo)}
\format{
  A data frame with 3 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{a factor with levels \code{P1} \code{P2} \code{P3}}
    \item{\code{STATION_NAME}}{a factor with levels \code{P1} \code{P2} \code{P3}}
    \item{\code{BASIN_NAME}}{a factor with levels \code{Basin}}
    \item{\code{Type}}{a factor with levels \code{Rain gauge}}
    \item{\code{X}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{Elevation}}{a numeric vector}
  }
}
\examples{
data(gisInfo)
## maybe str(gisInfo) ; plot(gisInfo) ...
}
\keyword{datasets}
