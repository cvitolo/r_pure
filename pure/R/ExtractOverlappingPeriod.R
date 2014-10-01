#' This function trims P, E and Q time series to the overlapping period.
#' It is possible to ignore Q by setting ignoreQ=TRUE
#'
#' @param tsList ist of time series (e.g. zoo objects)
#'
#' @return list of time series with the same temporal window
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # ExtractOverlappingPeriod(tsList)
#'


ExtractOverlappingPeriod <- function(tsList){

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

  x <- do.call(merge,newList)

  names(x) <- names(tsList)

  return(x)

}
