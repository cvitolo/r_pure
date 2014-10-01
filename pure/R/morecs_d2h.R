#' This function imports evapotranspiration data from MORECS (daily) and interpolate to hourly resolution
#'
#' @author Wouter Buytaert (2012/01/31), modified by Claudia Vitolo (2013/11/27)
#'
#' @param filepath string representing the file location, e.g. "~/Dropbox/Data/pontbren/ET/DAILY L Vyrnwy.CSV"
#'
#' @return zoo object containing the hourly time series of potential evapotranspiration
#'
#' @examples
#' #morecs_d2h( filepath="~/Dropbox/Data/pontbren/ET/DAILY L Vyrnwy.CSV" )
#'

morecs_d2h <- function(filepath){

  data <- read.csv(filepath)
  time <- as.POSIXct(strptime(paste(data[,1], data[,2], data[,3]),
                              format="%d %m %Y", tz="GMT"))

  ET <- data[,45]   # header = X.PE_GRS_M.

  ## disaggregate with sinusoidal function:

  ## sample a quadratic sinusoidal function:
  diurnal <- (0.5*sin((c(1:24)/24-0.25)*pi*2)+0.5)^2

  ## normalise:
  diurnal <- diurnal/sum(diurnal)

  ## matrix multiplication:
  EThourly <- as.matrix(ET) %*% diurnal

  time2 <- seq(time[1], time[length(time)] + 23*60*60, by="1 hour")

  EThourly <- zoo(as.vector(t(EThourly)), order.by=time2)

  return(EThourly)

}
