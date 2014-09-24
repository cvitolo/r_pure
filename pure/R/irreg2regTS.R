#' Transforms irregular time series to regular ones
#'
#' @param myTS can only be a single time series
#' @param deltim final time step
#'
#' @details Sometimes recording time in datasets can be shifted by few minutes.
#' This function forces the minutes to be multiple of a certain "multiplier".
#' It works with a maximum of 2 levels of nested time series. Something like:
#' list("P"=list(P1,P2,P3), "E"=E, "Q"=Q)
#'
#' @return time series with regular recording time
#'
#' @examples
#' # irr2reg( myTS, deltim )
#'

irr2reg <- function(myTS, deltim){

  multiplier=24*60*deltim
  newTS <- as.xts(aggregate(myTS, align.time(index(myTS), multiplier*60)))

  return(newTS)

}

#' Transforms irregular time series to regular ones
#'
#' @param CatchmentName string defining the name of the catchment
#' @param DataList number of minutes for the resulting regular time series. E.g. 15 (min)
#' @param deltim final time step
#'
#' @details Sometimes recording time in datasets can be shifted by few minutes.
#' This function forces the minutes to be multiple of a certain "multiplier".
#' It works only with PURE data with a maximum of 2 levels of nested time series. Something like: list("P"=list(P1,P2,P3), "E"=E, "Q"=Q)
#'
#' @return time series with regular recording time
#'
#' @examples
#' # irreg2regTS( CatchmentName, DataList, deltim )
#'

irreg2regTS <- function(CatchmentName, DataList, deltim){

  multiplier=24*60*deltim

  if (CatchmentName == "Pontbren") {

    print("Correcting P",quote=FALSE)
    P1 <- list()
    for ( raingauges in 1:length(as.list(DataList$P)) ) {
      P1[[raingauges]] <- as.xts(aggregate(DataList$P[[raingauges]],
                                           align.time(index(DataList$P[[raingauges]]),
                                                      multiplier*60)))
    }

    print("Correcting E",quote=FALSE)
    E1 <- as.xts(aggregate(DataList$E,
                           align.time(index(DataList$E),
                                      multiplier*60)))

    print("Correcting Q",quote=FALSE)
    Q1 <- as.xts(aggregate(DataList$Q,
                           align.time(index(DataList$Q),
                                      multiplier*60)))

  }
  if (CatchmentName == "Eden") {
    print(paste("No correction necessary for ", CatchmentName, " catchment."),quote=FALSE)
    P1 <- DataList$P
    E1 <- DataList$E
    Q1 <- DataList$Q
  }

  return(list("P"=P1,"E"=E1,"Q"=Q1))

}
