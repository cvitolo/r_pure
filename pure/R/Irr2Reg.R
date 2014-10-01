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
  myUnits <- attributes(median( diff( index(dataset) )))$units

  # calculate the multiplier in seconds
  if ( myUnits == "secs"  ) multiplier <- timestep
  if ( myUnits == "mins"  ) multiplier <- timestep*60
  if ( myUnits == "hours" ) multiplier <- timestep*60*60
  if ( myUnits == "days"  ) multiplier <- timestep*60*60*24

  newTS <- as.xts(aggregate(dataset, align.time(index(dataset), multiplier)))

  return(newTS)

}
