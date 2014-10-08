#' Find coarser temporal resolution amongst multiple time series
#'
#' @details This function finds the coarser temporal resolution amongst a list of time series.
#'
#' @param myList list containing multiple time series (e.g. zoo object)
#'
#' @return temporal resolution in seconds.
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # overlappingTS(myList)
#'

CommonTemporalResolution <- function(myList){

  # require(zoo)

  multiplier <- c()

  for (dataset in myList){

    # detect the most probable time step
    timestep <- median( as.numeric(diff( index(dataset) )) )
    myUnits <- attributes(median( diff( index(dataset) )))$units

    # calculate the multiplier in seconds
    if (myUnits == "secs")  multiplier <- append(multiplier, timestep)
    if (myUnits == "mins")  multiplier <- append(multiplier, timestep*60)
    if (myUnits == "hours") multiplier <- append(multiplier, timestep*60*60)
    if (myUnits == "days")  multiplier <- append(multiplier, timestep*60*60*24)

  }

  maxMultiplier <- max(multiplier)

  newList <- sapply(names(myList), function(x) NULL)
  counter <- 0

  for (dataset in myList){

    counter <- counter + 1
    m <- maxMultiplier/multiplier[counter]

     newList[[counter]] <- aggregate(dataset,
                                     align.time(index(dataset), maxMultiplier),
                                     FUN = mean)*m

  }

  names(newList) <- names(myList)

  return( list("maxMultiplier" = maxMultiplier, "aggregatedList" = newList) )

}
