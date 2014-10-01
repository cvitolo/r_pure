#' This function trims P, E and Q time series to the overlapping period.
#' It is possible to ignore Q by setting ignoreQ=TRUE
#'
#' @param P0 precipitation time series (e.g. zoo object)
#' @param E0 potential evapotranspiration time series (e.g. zoo object)
#' @param Q0 streamflow discharge time series (e.g. zoo object)
#' @param ignoreQ boolean by default set to FALSE. Used to the algorithm to ignore Q.
#'
#' @return zoo object in which P, E and Q are merged.
#'
#' @author Claudia Vitolo
#'
#' @examples
#' #ExtractOverlappingPeriod(P0,E0,Q0,ignoreQ=FALSE)
#'


ExtractOverlappingPeriod <- function(P0,E0,Q0,ignoreQ=FALSE){

  if (ignoreQ==TRUE) {
    # find starting time
    start0 <- max(index(P0)[1],index(E0)[1])
    # find end time
    end0 <- min(index(P0)[length(P0)],index(E0)[length(E0)])
  }else{
    # find starting time
    start0 <- max(index(P0)[1],index(E0)[1],index(Q0)[1])
    # find end time
    end0 <- min(index(P0)[length(P0)],index(E0)[length(E0)],index(Q0)[length(Q0)])
  }

  P <- window(P0, start=start0,end=end0)
  E <- window(E0, start=start0,end=end0)
  Q <- window(Q0, start=start0,end=end0)

  finalDATA <- zoo( merge(P,E,Q) )
  names(finalDATA) <- c("P","E","Q")

  return(finalDATA)
}
