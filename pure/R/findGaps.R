#' Finds gaps in time series
#'
#' @param myTS is a list of nested time series.
#' @param deltim final time step
#'
#' @return a list of fullranges, times and intervals of occurring gaps
#'
#' @examples
#' # findGaps( myTS, deltim )
#'

findGaps <- function(myTS,deltim) {

  multiplier <- 24*60*deltim

  values <- myTS$E
  datetime <- index(values)
  # Get the range of dates covered
  fullrangeE <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
  # Get the dates in fullrangeE that are not in DF$V2
  gapsE <- fullrangeE[!fullrangeE %in% datetime]
  x <- data.frame(datetime, values)
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
      values <- myTS$P[[raingauges]]
      datetime <- index(values)
      # Get the range of dates covered
      fullrangeP[[raingauges]] <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
      # Get the dates in fullrangeP that are not in DF$V2
      gapsP[[raingauges]] <- fullrangeP[[raingauges]][!fullrangeP[[raingauges]] %in% datetime]
      x <- data.frame(datetime, values)
      dtP[[raingauges]] <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
      if ( length(gapsP[[raingauges]])==0 ) {
        print(paste("No gaps in the time series P",raingauges,sep=""))
      }else{
        print(paste("Gaps in the time series P",raingauges," have been stored in the object gapsP[[",raingauges,"]]",sep=""))
      }
    }
  }else{
    values <- myTS$P
    datetime <- index(values)
    # Get the range of dates covered
    fullrangeP <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
    # Get the dates in fullrangeP that are not in DF$V2
    gapsP <- fullrangeP[!fullrangeP %in% datetime]
    x <- data.frame(datetime, values)
    dtP <- data.frame( start= x[c(1, diff(x$datetime))>1, ]$datetime, end=x[diff(x$datetime)>1, ]$datetime)
    if ( length(gapsP)==0 ) {
      print("No gaps in the time series P")
    }else{
      print("Gaps in the time series P have been stored in the object gapsP")
    }
  }

  values <- myTS$Q
  datetime <- index(values)
  # Get the range of dates covered
  fullrangeQ <- seq(min(datetime), max(datetime), by = 60*multiplier) # by seconds
  # Get the dates in fullrangeQ that are not in DF$V2
  gapsQ <- fullrangeQ[!fullrangeQ %in% datetime]
  x <- data.frame(datetime, values)
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
