\name{ScanTS}
\alias{ScanTS}
\title{This function scans the time series given in input and checks for:
- un-even recording times
- negative values
- missing values}
\usage{
  ScanTS(OriginalDataset, returnGapsInfo = TRUE,
    returnTimeInfo = TRUE, returnNegInfo = TRUE,
    verbose = TRUE)
}
\arguments{
  \item{OriginalDataset}{Time series, it can be zoo, xts,
  etc}

  \item{returnGapsInfo}{bolean, default value is TRUE}

  \item{returnTimeInfo}{bolean, default value is TRUE}

  \item{returnNegInfo}{bolean, default value is TRUE}

  \item{verbose}{bolean, default value is TRUE}
}
\value{
  print useful information
}
\description{
  This function scans the time series given in input and
  checks for: - un-even recording times - negative values -
  missing values
}
\examples{
#ScanTS( dataset=Q[[8]] , returnGapsInfo=TRUE, returnTimeInfo=TRUE, returnNegInfo=TRUE , verbose = TRUE)
}

