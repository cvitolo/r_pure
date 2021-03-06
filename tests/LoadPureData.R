#' Load data for any of the PURE case studies (13 + 2).
#'
#' @param CatchmentName this is a one-word string containing the name of the
#'        catchment. Accepted values are: "Eden" or "Pontbren"
#'
#' @param SubcatchmentName this is an integer number representing the name of
#'        the subcatchment. Predefined values are:
#'        For Eden:      1 (Eden at BlindBeck),   2 (Eden at Dacre),
#'        For Pontbren:  1 (Pontbren at Site1),   2 (Pontbren at Site2),
#'                       3 (Pontbren at Site3),   4 (Pontbren at Site4),
#'                       5 (Pontbren at Site5),   6 (Pontbren at Site6),
#'                       7 (Pontbren at Site7),   8 (Pontbren at Site8),
#'                       9 (Pontbren at Site9),  10 (Pontbren at Site10),
#'                      11 (Pontbren at Site11), 12 (Pontbren at Site12),
#'                      13 (Pontbren at Site13)
#'
#' @param datafolder path to the folder containing the data
#'
#' @return a list of 4 objects:
#'         $DATA a data.frame with 4 columns containing observations:
#'         date&time (Dates), discharge (Q), precipitation (P),
#'         potential evapo-transpiration (E)
#'         $deltimP observed precipitation time step
#'         $deltimE observed potential evapotranspiration time step
#'         $deltimQ observed streamflow discharge time step
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # LoadPureData(CatchmentName="Pontbren",SubcatchmentName=9, datafolder="/home/data/")
#'

LoadPureData <- function(CatchmentName,SubcatchmentName, datafolder){

  # LOAD DATA
  if (CatchmentName == "Pontbren") {

    load(paste(datafolder,"pontbren.rda",sep=""))
    pontbren <- pontbren

    E <- pontbren$Epontbren # 1 hour time step # as.numeric(median(diff(index(E))))
    P0 <- pontbren$Ppontbren # 10 min time step # E.g. as.numeric(median(diff(index(P[[1]]))))
    Q0 <- pontbren$Qpontbren # 15 min time step # E.g. as.numeric(median(diff(index(Q[[1]]))))

    if (SubcatchmentName == 1) {
      P <- list("P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 1311 # Total Length in m
      S <- 0.012 # Average slope in %
      A <- 1
      Area <- 0.630
    }

    if (SubcatchmentName == 2) {
      P <- list("P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 1700 # Total Length in m
      S <- 0.012 # Average slope in %
      A <- 1
      Area <- 1.198
    }

    if (SubcatchmentName == 3) {
      P <- list("P1"=P0[[1]],"P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 408 # Total Length in m
      S <- 0.011 # Average slope in %
      area <- 0.263
      A <- c(0.093/area,0.17/area)
      Area <- 0.263
    }

    if (SubcatchmentName == 4) {
      P <- list("P1"=P0[[1]],"P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 3822 # Total Length in m
      S <- 0.015 # Average slope in %
      area <- 2
      A <- c(0.632/area,1.368/area)
      Area <- 2
    }

    if (SubcatchmentName == 5) {
      P <- list("P1"=P0[[1]],"P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 3925 # Total Length in m
      S <- 0.016 # Average slope in %
      area <- 2.4
      A <- c(0.873/area,1.527/area)
      Area <- 2.4
    }

    if (SubcatchmentName == 6) {
      P <- list("P1"=P0[[1]],"P4"=P0[[4]],"P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 5603 # Total Length in m
      S <- 0.017 # Average slope in %
      area <- 3.171
      A <- c(1.464/area,0.178/area,1.529/area)
      Area <- 3.171
    }

    if (SubcatchmentName == 7) {
      P <- list("P1"=P0[[1]],"P2"=P0[[2]],"P3"=P0[[3]],"P4"=P0[[4]],"P5"=P0[[5]],"P6"=P0[[6]])
      Q <- Q0[[SubcatchmentName]]
      L <- 6738 # Total Length in m
      S <- 0.021 # Average slope in %
      area <- 5.718
      A <- c(2.454/area,0.260/area,0.119/area,0.323/area,2.449/area,0.113/area)
      Area <- 5.718
    }

    if (SubcatchmentName == 8) {
      P <- list("P3"=P0[[3]])
      Q <- Q0[[SubcatchmentName]]
      L <- 2238 # Total Length in m
      S <- 0.002 # Average slope in %
      A <- 1
      Area <- 1.288
    }

    if (SubcatchmentName == 9) {
      P <- list("P2"=P0[[2]],"P3"=P0[[3]],"P5"=P0[[5]])
      Q <- Q0[[SubcatchmentName]]
      L <- 8097 # Total Length in m
      S <- 0.015 # Average slope in %
      area <- 4.083
      A <- c(1.347/area,2.722/area,0.014/area)
      Area <- 4.083
    }

    if (SubcatchmentName == 10) {
        P <- list("P1"=P0[[1]],"P2"=P0[[2]],"P3"=P0[[3]],"P4"=P0[[4]],"P5"=P0[[5]])
        Q <- NA # Q0[[SubcatchmentName]]
        print("Discharge data not available at this location.")
        L <- 11410 # Total Length in m
        S <- 0.019 # Average slope in %
        area <- 12.507
        A <- c(2.551/area,3.232/area,2.872/area,0.349/area,2.456/area,1.047/area)
        Area <- 12.507
    }

    if (SubcatchmentName == 11) {
      P <- list("P6"=P0[[6]])
      Q <- NA # Q0[[SubcatchmentName]]
      print("Discharge data not available at this location.")
      L <- 1695 # Total Length in m
      S <- 0.019 # Average slope in %
      A <- 1
      Area <- 0.765
    }

    if (SubcatchmentName == 12) {
      P <- list("P6"=P0[[6]])
      Q <- NA # Q0[[SubcatchmentName]]
      print("Discharge data not available at this location.")
      L <- 1755 # Total Length in m
      S <- 0.020 # Average slope in %
      A <- 1
      Area <- 1.048
    }

    if (SubcatchmentName == 13) {
      P <- list("P6"=P0[[6]])
      Q <- Q0[[SubcatchmentName]]
      L <- 2830 # Total Length in m
      S <- 0.019 # Average slope in %
      A <- 1
      Area <- 1.290
    }

    deltimP <- 1/24/6
    deltimE <- 1/24
    deltimQ <- 1/24/4

    #DTM <- readGDAL(paste(datafolder,"Wales_Terra50_merged.tif",sep=""))
    #DTM <- readGDAL(paste(datafolder,"PontbrenTerra50.tif",sep=""))
    #Soil <- readGDAL(paste(datafolder,"PontbrenSoilraster.tif",sep=""))
    #LandUse <- readGDAL(paste(datafolder,"Pontbrenlcm2000grd.tif",sep=""))
    #LandUse <- raster(paste(datafolder,"Pontbrenlcm2000grd.tif",sep=""))

    DTM <- raster(paste(datafolder,"PontbrenTerra50.tif",sep=""))
    Soil <- raster(paste(datafolder,"PontbrenSoilraster.tif",sep=""))
    LandUse <- raster(paste(datafolder,"LCM2007_pontbren.tif",sep=""))
  }

  if (CatchmentName == "Eden") {

    load(paste(datafolder,"eden.rda",sep=""))
    eden <- eden

    if (SubcatchmentName == 1) {  #BlindBeck
      # It uses Blindbeck dataset from Sep 2005 to Oct 2005
      DATA0 <- eden$blinbeckdata
      P <- zoo(DATA0[,"RAIN"],order.by=DATA0[,1])
      E <- zoo(DATA0[,"PET"],order.by=DATA0[,1])
      Q <- zoo(DATA0[,"FLOW"],order.by=DATA0[,1])
      L <- 5961 # Total Length in m
      S <- 0.032  # Average slope in %
      A <- 1
      Area <- 10.332
      DTM <- raster(paste(datafolder,"BlindBeck_Terra50_merged.tif",sep=""))
      #Soil <- raster(paste(datafolder,"BlindBeck_Soils.tif",sep=""))
      LandUse <- raster(paste(datafolder,"LCM2007_blindbeck.tif",sep=""))
    }

    if (SubcatchmentName == 2) { # Dacre
      # It uses Dacre dataset from 01/12/2004 to 16/01/2005
      DATA0 <- eden$dacredata
      P <- zoo(DATA0[,"RAIN"],order.by=DATA0[,1])
      E <- zoo(DATA0[,"PET"],order.by=DATA0[,1])
      Q <- zoo(DATA0[,"FLOW"],order.by=DATA0[,1])
      L <- 4578 # Total Length in m
      S <- 0.017 # Average slope in %
      A <- 1
      Area <- 9.880
      DTM <- raster(paste(datafolder,"Dacre_Terra50_merged.tif",sep=""))
      #Soil <- raster(paste(datafolder,"Dacre_Soils.tif",sep=""))
      LandUse <- raster(paste(datafolder,"LCM2007_dacre.tif",sep=""))
    }

    deltimP <- deltimE <- deltimQ <- 1/24/4
  }

  if (CatchmentName != "Pontbren" & CatchmentName != "Eden") {
    print("Invalid name of catchment! See documentation.")
  }else{

    ID <- NULL

    # Projection for the GIS files
    # British National Grid
    # proj4string(DTM) <- CRS("+init=epsg:27700")
    GISinfo <- read.csv( paste(datafolder,"pure_rain_stations.csv",sep="") )

    # do not use ~ for home folder, rgdal does not like it!
    subcatchments <- readOGR(dsn=datafolder,layer="subcatchments")

    catchment <- subset(subcatchments,
                        ID==paste(ifelse(CatchmentName=="Pontbren","P","E"),
                                  SubcatchmentName,sep=""))

    DATA <- list("P"=P,
                 "E"=E,
                 "Q"=Q,
                 "deltimP"=deltimP ,
                 "deltimE"=deltimE,
                 "deltimQ"=deltimQ,
                 "L"=L,
                 "S"=S,
                 "A"=A,
                 "Area"=Area,
                 "DTM"=DTM,
                 "Soil"=Soil,
                 "LandUse"=LandUse,
                 "subcatchments"=subcatchments,
                 "catchment" = catchment,
                 "GISinfo"=GISinfo)
  }

  return(DATA)

}
