#' Function to correct negative values.
#'
#' @param regTSList List of datasets containing P, E and Q (with corrected recording times)
#'
#' @return a dataset with non negative values
#'
#' @details At the moment the function simply substitutes negative values with NAs.
#'
#' @examples
#' # correctNeg( regTSList )
#'

correctNeg <- function(regTSList){

  P1 <- regTSList$P

  for ( raingauges in 1:length(as.list(P1)) ) {
    if ( any(P1[[raingauges]]<0) ) {
      print(paste("Correction is necessary for P time series"),quote=FALSE)
      P1[which(P1<0)] <- NA
    }else{
      print(paste("No correction necessary for P",raingauges," time series",sep=""),quote=FALSE)
    }
  }

  print(paste("No correction necessary for E time series"),quote=FALSE)

  Q1 <- regTSList$Q
  if ( any(Q1<0) ) {
    print(paste("Correction is necessary for Q time series"),quote=FALSE)
    Q1[which(Q1<0)] <- NA
  }else{
    print(paste("No correction necessary for Q time series"),quote=FALSE)
  }

  return(list("P"=P1,"E"=regTSList$E,"Q"=Q1))

}
