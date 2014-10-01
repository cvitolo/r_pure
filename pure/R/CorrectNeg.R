#' Correct negative values.
#'
#' @param tseries time series
#'
#' @return a time series with non negative values
#'
#' @details At the moment the function simply substitutes negative values with NAs.
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # CorrectNeg( tseries )
#'

CorrectNeg <- function( tseries ){

  if ( any( tseries < 0 ) ) {

      message("I found negative values and changed them to NA.")
      message("It is now necessary to infill the generated missing values.")

      negativeElements <- which(tseries<0)
      tseries[negativeElements] <- NA

    }else{

      message("No negative values found.")

    }

  return(tseries)

}
