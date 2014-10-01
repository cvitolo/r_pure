#' Calculate Potential Evapotranspiration using Penman Monteith formula
#'
#' @param stationElevation elevation of the recording station (in metres)
#' @param TD time series of dry bulb temperature
#' @param TW time series of wet bulb temperature
#' @param NR time series of net radiation
#' @param WS time series of wind speed
#'
#' @return ET0 = ET for grass (K = 1)
#'
#' @references FAO - hourly method: http://www.fao.org/docrep/x0490e/x0490e08.htm#TopOfPage
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # data(exampleTimeSeries)
#' # pet(stationElevation=0,TD,TW,NR,WS)
#'

pet <- function(stationElevation=0,TD,TW,NR,WS){

  # Net Radiation usually needs transformation
  NR <- NR/(1000000/36000)

  TDtemp <- coredata(TD)
  TWtemp <- coredata(TW)
  NRtemp <- coredata(NR)
  WStemp <- coredata(WS)

  myHOUR <- as.POSIXlt(index(NR))$hour

  G <- rep(0,length(NRtemp))

  for (r in 1:length(NRtemp)){

    # print(paste(r,"out of",length(NRtemp)))

    if (myHOUR[r] >= 9 & myHOUR[r] <= 17) {
      G[r] <- 0.1 * NRtemp[r]
    }
    if ( (myHOUR[r] < 9 & myHOUR[r] >= 0) | (myHOUR[r] > 17 & myHOUR[r] <= 24) ) {
      G[r] <- 0.5 * NRtemp[r]
    }
  }

  delta <- as.vector(4098*(0.6108*exp(17.27*TDtemp/(TDtemp+237.3)))/(TDtemp+237.3)^2)

  satvapourpressure <- 0.6108*exp(17.27*TWtemp/(TWtemp+237.3))

  # Psychometric constant (gam) for different altitudes (z)
  # (Source: FAO Irrigation and Drainage Paper No.56, p.214)
  elevation <- seq(from=0, to=4000, by=100) # in metres
  gam <- c(0.067,0.067,0.066,0.065,0.064,0.064,0.063,0.062,0.061,0.061,
           0.060,0.059,0.058,0.058,0.057,0.056,0.056,0.055,0.054,0.054,
           0.053,0.052,0.052,0.051,0.051,0.050,0.049,0.049,0.048,0.047,
           0.047,0.046,0.046,0.045,0.045,0.044,0.043,0.043,0.042,0.042,0.041)
  known <- data.frame(elevation,gam)
  #plot (known$elevation, known$gam, type="o")
  model.lm <- lm(gam ~ elevation, data = known)
  # Use predict to estimate the values for stationElevation.
  # Note that predict expects a data.frame and the col
  # names need to match
  newY <- predict(model.lm, newdata = data.frame(elevation = stationElevation))
  #Add the predicted points to the original plot
  #points(stationElevation, newY, col = "red")

  # The psychometric constant for Severn is newY
  actualvapourpressure <- satvapourpressure-newY*(TDtemp-TWtemp)

  Etemp <- as.vector(0.408*delta*(NRtemp-G)+newY*900/(TDtemp+273)*WStemp*(satvapourpressure-actualvapourpressure)/(delta+newY*(1+0.34*WStemp)))

  E <- zoo(Etemp,order.by=index(NR))
  #plot(E);max(E,na.rm=TRUE)

  return(round(E,2))

}
