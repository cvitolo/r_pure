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
#' # Irr2Reg( dataset, deltim )
#'

Irr2Reg <- function(dataset){

  # detect the most probable time step
  timestep <- median( as.numeric(diff( index(dataset) )) )
  myUnits <- attributes(median( diff( index(dataset) )))$units

  if ( myUnits == "mins"  ) deltim <- timestep
  if ( myUnits == "hours" ) deltim <- timestep*60
  if ( myUnits == "days"  )  deltim <- timestep*60*24

  multiplier <- 60*deltim # seconds
  newTS <- as.xts(aggregate(dataset, align.time(index(dataset), multiplier)))

  return(newTS)

}
