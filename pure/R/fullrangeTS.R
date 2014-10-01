#' Create a time series without gaps in recording
#'
#' @param myTS is a list of nested time series.
#' @param myRanges range obtained from function \code{findGaps}
#'
#' @return time series with no gaps in recording
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # fullrangeTS( myTS, myRanges )
#'

fullrangeTS <- function(myTS, myRanges){

  values <- myTS$E
  range <- myRanges$fullrangeE
  dummy <- zoo(rep(NA,length(range)),order.by=range)
  temp <- merge(values, dummy, all = TRUE, fill = NA)[,"values"]
  mergedE <- zoo(as.numeric(temp),order.by=index(temp))

  if ( typeof(myTS$P)=="list" ){

    mergedP <- list()

    for ( raingauges in 1:length(as.list(myTS$P)) ) {

      values <- myTS$P[[raingauges]]
      range <- myRanges$fullrangeP[[raingauges]]
      dummy <- zoo(rep(NA,length(range)),order.by=range)
      temp <- merge(values, dummy, all = TRUE, fill = NA)[,"values"]
      mergedP[[raingauges]] <- zoo(as.numeric(temp),order.by=index(temp))

    }

  }else{

    values <- myTS$P
    range <- myRanges$fullrangeP
    dummy <- zoo(rep(NA,length(range)),order.by=range)
    temp <- merge(values, dummy, all = TRUE, fill = NA)[,"values"]
    mergedP <- zoo(as.numeric(temp),order.by=index(temp))

  }

  values <- myTS$Q
  range <- myRanges$fullrangeQ
  dummy <- zoo(rep(NA,length(range)),order.by=range)
  temp <- merge(values, dummy, all = TRUE, fill = NA)[,"values"]
  mergedQ <- zoo(as.numeric(temp),order.by=index(temp))

  merged <- list("P"=mergedP,"E"=mergedE,"Q"=mergedQ)

  return(merged)

}
