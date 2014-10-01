#' Calculate the areal average based on a list of time series (e.g. from rain gauges)
#'
#' @param tsList list containing time series
#' @param areas percentage of total area covered by the Voronoi polygons related to each time series
#' @param gisInfo data.frame containing 6 variables: ID, STATION_NAME, BASIN_NAME,Type, X,Y, Elevation
#' @param catchmentShape SpatialPolygonsDataFrame containing the catchment boundary
#' @param interpolationMethod default value is "Thiessen". Available algorithms are: "Mean", "Thiessen",  "IDW" (Inverse Distance Weighted), "OK" (Ordinary Kriging)
#'
#' @author Claudia Vitolo
#'
#' @return areal averaged time series
#'

ArealAveraging <- function(tsList, areas=NULL,
                           gisInfo=NULL, catchmentShape=NULL,
                           interpolationMethod = "Thiessen"){

  timeIndex <- tsList[,1]
  mat <- tsList[,2:dim(tsList)[2]]

  if ( interpolationMethod =="Mean" ) {

    arealAveragedTS <- apply(mat,1,mean)

  }

  if ( interpolationMethod =="Thiessen" ) {

    # Specify the areas to perform areal average
    # (areas = percentages of total area)
    arealAveragedTS <- apply(mat, 1, weighted.mean, w = areas, na.rm = TRUE)

  }

  # INVERSE DISTANCE
  if ( interpolationMethod =="IDW" ) {

    # gstat does not like missing data, subset original data:
    # d <- na.omit(d)
    # convert simple data frame into a spatial data frame object:
    coordinates(gisInfo) <- ~ X+Y

    arealAveragedTS <- hydrokrige(x.ts = mat,
                                  x.gis = data.frame(gisInfo),
                                  X ="X", Y = "Y", sname = "ID",
                                  type = "block",
                                  dates = timeIndex,
                                  subcatchments = catchmentShape,
                                  plot = FALSE,
                                  verbose = FALSE)
  }

  # ORDINARY KRIGING
  if ( interpolationMethod =="OK" ) {

    # gstat does not like missing data, subset original data:
    # d <- na.omit(d)
    # convert simple data frame into a spatial data frame object:
    coordinates(gisInfo) <- ~ X+Y

    arealAveragedTS <- hydrokrige(x.ts = mat,
                                  x.gis = data.frame(gisInfo),
                                  X = "X", Y = "Y", sname = "ID",
                                  type = "block",
                                  dates = timeIndex,
                                  formula = value~1,
                                  subcatchments = catchmentShape,
                                  plot = FALSE,
                                  verbose = FALSE)

  }

  return( zoo(arealAveragedTS, order.by=timeIndex) )

}
