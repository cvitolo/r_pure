#' Identify Rainfall-Runoff events
#'
#' @param data data.frame of precipitation (P) and discharge (Q) time series
#' @param dt data time step
#'
#' @author Claudia Vitolo
#'
#' @return areal averaged time series
#'

EventIdentification <- function(data,dt){

  # require(zoo)
  # require(EcoHydRology)

  # STEP 1
  # Separate baseflow from runoff
  bfs <- zoo(BaseflowSeparation(coredata(data$Q))[,"bt"], order.by=index(data))
  ## You can check out how this looks with the hydrograph function:
  # hydrograph(input=data.frame("datetime"=index(data),"P"=coredata(data$P),
  #                             "Q"=coredata(data$Q),"Qbase"=coredata(bfs[,1])))
  # Calculate the median of baseflow to use as threshold
  threshold <- median(bfs)
  # Use quick flow only (runoff - baseflow)
  qfl <- zoo(BaseflowSeparation(coredata(data$Q))[,"qft"], order.by=index(data))

  # STEP 2
  # Identify discrete RAINFALL events from time series using hydromad::eventseq.
  evp <- eventseq(data$P,        # precipitation time series
                  thresh = 0.1,  # threshold value
                  inthresh = 0,  # second threshold to define when events stop
                  indur = 6)     # P must remain below inthresh for indur time
                                 # steps in order to terminate an event
  ## You can check out the structure of the evp object
  # str(evp)

  # STEP 3
  # Calculate rainfall centroid
  numberOfEvents <- dim(eventinfo(data$P,evp))[1]
  for (n in 1:numberOfEvents) {
    start <- eventinfo(data$P,evp)$Time[n]
    end   <- eventinfo(data$P,evp)$Time[n] + eventinfo(data$P,evp)$Duration[n]*x
    sumP  <- /sum(coredata(window(data$P,start=start,end=end)))
  }

  # Extend the events to include RUNOFF recession limb.
  events <- data.frame("Start"     = eventinfo(data$P,evp)$Time,
                       "End"       = rep(NA,numberOfEvents),
                       "CentroidP" = rep(NA,numberOfEvents),
                       "CentroidQ" = rep(NA,numberOfEvents),
                       "")

  # STEP 3
  # Calculate centroids of rainfall and runoff
  for (event in 1:dim(eventinfo(data$P,evp))[1]){

    eventinfo(data$P,evp)$Time

  }
  ## You can check out how this looks with the hydrograph function:
  # hydrograph(input=data.frame("datetime"=index(data), "P"=coredata(data$P),
  #                             "Q"=coredata(data$Q), "Qbase"= bfs[,1]))



}
