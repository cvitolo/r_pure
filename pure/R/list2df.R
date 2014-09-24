#' Convert a list of time series (from rain gauges) to a data frame
#'
#' @param P list of objects containing rainfall gauges time series
#'
#' @return data frame object
#'
#' @details P should contain objects with same length, start and end date-time. Use function overlappingTS. See package hydroTSM for details on the x.ts and x.gis objects.
#'
#' @examples
#' # list2df(P)
#'

list2df <- function(P){

  Pdf <- data.frame(matrix(unlist(P), ncol=length(P)))
  colnames(Pdf) <- names(P)
  PdfwithDateTime <- cbind("date-time"=index(P[[1]]),Pdf)

  return(PdfwithDateTime)

}
