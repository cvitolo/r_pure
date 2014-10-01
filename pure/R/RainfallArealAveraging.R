#' Function to calculate the rainfall areal average based on a list of time series (from rain gauges)
#'
#' @param CatchmentName string identifying the name of the catchment
#' @param Plist list containing rainfall gauges time series
#' @param DataList list of objects from function LoadMyData
#' @param InterpolationMethod default value is "Thiessen". Available algorithms are: "AritmeticMean", "Thiessen",  "IDW" (Inverse Distance Weighted), "OK" (Ordinary Kriging)
#'
#' @author Claudia Vitolo
#'
#' @return areal averaged time series
#'

RainfallArealAveraging <- function(CatchmentName, Plist, DataList, InterpolationMethod = "Thiessen"){

  RainfallDF <- list2df(Plist)

  if ( length(Plist) == 1 ) {

    P0 <- Plist$P[[1]]

  }else{

    print(paste("Apply choosen INTERPOLATION method: ",InterpolationMethod, sep=""),quote=FALSE)

    l <- length(Plist) + 1

    if ( InterpolationMethod =="AritmeticMean" ) {

      P0 <- apply(RainfallDF[,2:l],1,mean)

    }

    if ( InterpolationMethod =="Thiessen" ) {

      # Specify the areas to perform areal average (Areas = percentage of total area)
      A <- rep(DataList$A[1],dim(RainfallDF)[1])
      for ( counter in 2:(l-1) ) {
        A <- cbind(A,rep(DataList$A[counter],dim(RainfallDF)[1]))
      }

      P0 <- c()
      for ( counter in 2:l ) {
        P1 <- cbind(P0,RainfallDF[,counter]*A[,counter-1])
        P0 <- apply(P1,1,sum)
      }

    }

    if ( InterpolationMethod =="IDW" | InterpolationMethod =="OK" | InterpolationMethod =="KED" ) {

      ## load data
      GISinfo <- DataList$GISinfo ## gstat does not like missing data, subset original data: d <- na.omit(d)

      ## extract values based on catchment name
      if (CatchmentName == "Pontbren") RaingaugesLocations <- GISinfo[1:6,]
      if (CatchmentName == "Eden") RaingaugesLocations <- GISinfo[7:8,]

      ## convert simple data frame into a spatial data frame object:
      coordinates(RaingaugesLocations) <- ~ X+Y
      SubcatchmentName <- NULL
      shape <- DataList$subcatchments[DataList$subcatchments@data$ID == paste(ifelse(CatchmentName=="Pontbren","P","E"),SubcatchmentName,sep="") ,]

      # INVERSE DISTANCE
      if ( InterpolationMethod =="IDW" ) {

        P0 <- rep(NA, dim(RainfallDF)[1])

        for (t in 1:dim(RainfallDF)[1]) {

        P0[t] <- as.numeric(hydroTSM::hydrokrige(x.ts = RainfallDF[1,],
                                                 x.gis = DataList$GISinfo,
                                                 X ="X", Y = "Y", sname = "ID",
                                                 bname = "STATION_NAME",
                                                 type = "block",
                                                 subcatchments = shape,
                                                 cell.size = 50,
                                                 ColorRamp = "Precipitation",
                                                 main = "IDW Precipitation",
                                                 dates = RainfallDF[t,"date-time"],
                                                 verbose = FALSE))[2]
        }
      }

      # ORDINARY KRIGING
      if ( InterpolationMethod =="OK" ) {

        P0 <- rep(NA, dim(RainfallDF)[1])

        for (t in 1:dim(RainfallDF)[1]) {

          P0[t] <- as.numeric(hydroTSM::hydrokrige(x.ts = RainfallDF[1,],
                                                   x.gis = DataList$GISinfo,
                                                   X = "X", Y = "Y",
                                                   sname = "ID",
                                                   bname = "STATION_NAME",
                                                   type = "block",
                                                   formula = value~1,
                                                   subcatchments = shape,
                                                   cell.size = 50,
                                                   ColorRamp = "Precipitation",
                                                   main = "OK Precipitation",
                                                   arrow.plot = TRUE,
                                                   arrow.offset = c(900000,4750000),
                                                   arrow.scale = 20000,
                                                   scalebar.plot = TRUE,
                                                   sb.offset = c(400000,4480000),
                                                   sb.scale = 100000,
                                                   dates = RainfallDF[t,"date-time"],
                                                   verbose = FALSE))[2]
        }
      }

    }

  }

  return( zoo(P0, order.by=RainfallDF[,1]) )

}
