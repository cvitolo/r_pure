#' Trim all the time series contained in a list to the common period.
#'
#' @param tsList list of time series (e.g. zoo objects)
#'
#' @return list of time series with the same temporal window
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # CommonRecordingTime(tsList)
#'

CommonRecordingTime <- function(tsList){

  maxIndex <- c()
  minIndex <- c()

  for(myTS in tsList) {

    maxIndex <- append(maxIndex, max(index(myTS)))
    minIndex <- append(minIndex, min(index(myTS)))

  }

  start0 <- max(minIndex)
  end0   <- min(maxIndex)

  newList <- list()

  counter <- 0

  for(myTS in tsList) {

    counter <- counter + 1

    newList[[counter]] <- window(myTS, start=start0,end=end0)

  }

  names(newList) <- names(tsList)

  #x <- do.call(merge,newList)

  #names(x) <- names(tsList)

  return(newList)

}
