#' Transforms irregular time series to regular ones
#'
#' @param dataset can only be a single time series
#'
#' @details Sometimes recording time in datasets can be shifted by few minutes.
#' This function forces the minutes to be multiple of a certain "multiplier".
#'
#' @return time series with regular recording time
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # Irr2Reg( dataset )
#'

Irr2Reg <- function(dataset){

  # detect the most probable time step
  timestep <- median( as.numeric(diff( index(dataset) )) )
  myUnits  <- attributes(median( diff( index(dataset) )))$units

  startT <- index(dataset)[1]
  endT   <- index(dataset)[length(dataset)]
  difT   <- difftime(endT,startT)

  # calculate the multiplier in seconds

  if ( myUnits == "secs"  ) multiplier <- timestep
  if ( myUnits == "mins"  ) multiplier <- timestep*60
  if ( myUnits == "hours" ) multiplier <- timestep*60*60
  if ( myUnits == "days"  ) multiplier <- timestep*60*60*24

  myMin <- 60 - timestep
  myDif <- difT[[1]]*24*60/timestep

  # newTS <- as.xts(aggregate(dataset, align.time(index(dataset), multiplier)))

  tempDates <- ISOdatetime(as.POSIXlt(startT)$year + 1900,
                           as.POSIXlt(startT)$mon + 1,
                           as.POSIXlt(startT)$mday,
                           as.POSIXlt(startT)$hour,
                           myMin, # as.POSIXlt(startT)$min
                           0,     # as.POSIXlt(startT)$sec
                           ) + seq(0:(myDif-1)) * multiplier

  newTS <- na.locf(merge(xts(,tempDates),dataset))[tempDates]

  return(newTS)

}
