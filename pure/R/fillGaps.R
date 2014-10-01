#' Infill missing values.
#'
#' @param myList List of datasets containing P, E and Q (with corrected recording times and no negative values for P and Q)
#'
#' @return a dataset with infilled missing values
#'
#' @details It uses the fuction na.approx (better not use na.spline, as it generates negative values)
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # fillGaps( myList )
#'

fillGaps <- function(myList){

  E <- na.approx(myList$E)

  if ( typeof(myList$P)=="list"  ){
    P <- list()
    for ( raingauges in 1:length(as.list(myList$P)) ) {
      P[[raingauges]] <- na.approx(myList$P[[raingauges]])
    }
  }else{
    P <- na.approx(myList$P)
  }

  Q <- na.approx(myList$Q)

  return(list("P"=P,"E"=E,"Q"=Q))

}
