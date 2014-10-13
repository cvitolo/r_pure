#' This function trims a list made of multiple datasets (e.g. multiple rain gauges datsets) to the overlapping period.
#'
#' @param myList list containing multiple time series (e.g. zoo object)
#'
#' @return list of zoo objects in which each time series is defined for the same overlapping period.
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # OverlappingTS(myList)
#'

OverlappingTS <- function(myList){

  if ( typeof(myList$P)=="list"  ){
    P <- as.list(zoo(do.call("merge", Map(as.zoo, c(myList$P, all = FALSE)))))
    Psample <- zoo(P[[1]])
  }else{
    P <- Psample <- zoo(myList$P)
  }

  temp <- setNames(myList, c("P", "E", "Q"))
  List1 <- list("P"=Psample,"E"=zoo(myList$E),"Q"=zoo(myList$Q))
  List2 <- as.list(zoo(do.call("merge", Map(as.zoo, c(List1, all = FALSE)))))

  List2$P <- P

  return(List2)

#
#
#   # find starting time
#   start0 <- index(myList[[1]][1])
#   for ( i in 2:length(myList) ) {
#     if ( start0 < index(myList[[i]][1]) ) {
#       start0 <- index(myList[[i]][1])
#     }
#   }
#
#   # find end time
#   end0 <- index(myList[[1]][length(myList[[1]])])
#   for ( i in 2:length(myList) ) {
#     if ( end0 > index(myList[[i]][length(myList[[i]])]) ) {
#       end0 <- index(myList[[i]][length(myList[[i]])])
#     }
#   }
#
#   newmyList <- list()
#   #extract only overlapping periods
#   for ( i in 1:length(myList) ) {
#     newmyList[[i]] <- window(myList[[i]], start=start0,end=end0)
#   }
#
#   return(newmyList)

}
