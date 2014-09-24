#' This function converts inputs "from" one unit "to" another.
#'
#' @param dataset Time series, it can be zoo, xts, etc
#' @param from string to specify the original unit, e.g. "m3/s" or "l/s"
#' @param to string to specify the destination unit, e.g. "mm/day" (default units to use with FUSE)
#' @param optionalInput real number which is used in case of special convertions (e.g. in case discharge is given in (from=) l/s and this needs to be converted to (to=)mm, then additional input is the area of the catchment)
#'
#' @return converted dataset
#'
#' @references For the "from" and "to" strings, check the Unidata's udunits reference: http://www.unidata.ucar.edu/software/udunits/
#'
#' @examples
#' #ConvertTS( dataset, from, to, optionalInput=NA )
#'

ConvertTS <- function(dataset, from, to, optionalInput=NA){

  if ( (from=="l/s" | from =="m3/s") & (to=="mm/day") ) {

    print( paste("Note: to convert from",from, "to", to, "the 'optionalInput' is recognized as an area with unit compatible with the 'to' parameter") )

    modifiedDataset <- udunits2::ud.convert(dataset, from, "mm3/day")/optionalInput

  }else{

    if (from=="mm/15min" & to=="mm/day") {
      modifiedDataset <- dataset*4*24
    }else{

      if (from=="mm/10min" & to=="mm/day") {
        modifiedDataset <- dataset*6*24
      }else{
        modifiedDataset <- udunits2::ud.convert(dataset, from, to)
      }

    }

  }

  return(modifiedDataset)

}
