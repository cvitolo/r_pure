#' Finds gaps in time series
#'
#' @param myTS is a list of nested time series.
#' @param deltim final time step
#'
#' @return a list of fullranges, times and intervals of occurring gaps
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # findGaps( myTS, deltim )
#'

findGaps <- function(myTS,deltim) {

  multiplier <- 24*60*deltim

  dataset <- myTS$E
  datetime <- index(dataset)
  # Get the range of dates covered
  fullrangeE <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
  # Get the dates in fullrangeE that are not in DF$V2
  gapsE <- fullrangeE[!fullrangeE %in% datetime]
  x <- data.frame(datetime, dataset)
  dtE <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
  if ( length(gapsE)==0 ) {
    print("No gaps in the time series E")
  }else{
    print("Gaps in the time series E have been stored in the object gapsE")
  }

  if ( typeof(myTS$P)=="list" ){
    gapsP <- list()
    fullrangeP <- list()
    dtP <- list()
    for ( raingauges in 1:length(as.list(myTS$P)) ) {
      dataset <- myTS$P[[raingauges]]
      datetime <- index(dataset)
      # Get the range of dates covered
      fullrangeP[[raingauges]] <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
      # Get the dates in fullrangeP that are not in DF$V2
      gapsP[[raingauges]] <- fullrangeP[[raingauges]][!fullrangeP[[raingauges]] %in% datetime]
      x <- data.frame(datetime, dataset)
      dtP[[raingauges]] <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
      if ( length(gapsP[[raingauges]])==0 ) {
        print(paste("No gaps in the time series P",raingauges,sep=""))
      }else{
        print(paste("Gaps in the time series P",raingauges," have been stored in the object gapsP[[",raingauges,"]]",sep=""))
      }
    }
  }else{
    dataset <- myTS$P
    datetime <- index(dataset)
    # Get the range of dates covered
    fullrangeP <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
    # Get the dates in fullrangeP that are not in DF$V2
    gapsP <- fullrangeP[!fullrangeP %in% datetime]
    x <- data.frame(datetime, dataset)
    dtP <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
    if ( length(gapsP)==0 ) {
      print("No gaps in the time series P")
    }else{
      print("Gaps in the time series P have been stored in the object gapsP")
    }
  }

  dataset <- myTS$Q
  datetime <- index(dataset)
  # Get the range of dates covered
  fullrangeQ <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
  # Get the dates in fullrangeQ that are not in DF$V2
  gapsQ <- fullrangeQ[!fullrangeQ %in% datetime]
  x <- data.frame(datetime, dataset)
  dtQ <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
  if ( length(gapsQ)==0 ) {
    print("No gaps in the time series Q")
  }else{
    print("Gaps in the time series Q have been stored in the object gapsQ")
  }

  fullranges <- list("fullrangeP"=fullrangeP,"fullrangeE"=fullrangeE,"fullrangeQ"=fullrangeQ)
  times <- list("gapsP"=gapsP,"gapsE"=gapsE,"gapsQ"=gapsQ)
  intervals <- list("dtP"=dtP,"dtE"=dtE,"dtQ"=dtQ)

  return(list("fullranges"=fullranges,"times"=times,"intervals"=intervals))

}
