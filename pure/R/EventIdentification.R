#' Identify Rainfall-Runoff events
#'
#' @param data data.frame of precipitation (P) and discharge (Q) time series
#'
#' @author Claudia Vitolo
#'
#' @return data.frame containing 3 columns: Tr (return period), P (max precipitation) and Q (max discharge)
#'

EventIdentification <- function(data){

  # require(zoo)
  # require(EcoHydRology)
  # data <- DATA[1:500,c("P","Q")]

  ### STEP 1 ###
  # Separate baseflow from runoff
  bfs <- zoo(BaseflowSeparation(coredata(data$Q))[,"bt"], order.by=index(data))
  ## You can check out how this looks with the hydrograph function:
  # hydrograph(input=data.frame("datetime"=index(data),"P"=coredata(data$P),
  #                             "Q"=coredata(data$Q),"Qbase"=coredata(bfs[,1])))

  # Use quick flow only (runoff - baseflow)
  qfl <- zoo(BaseflowSeparation(coredata(data$Q))[,"qft"], order.by=index(data))

  # set up a warmup period
  percentageWarmUp <- 10
  warmup <- round(dim(data)[1]/percentageWarmUp,0)
  pperiod <- (warmup+1):dim(data)[1]

  bfs <- bfs[pperiod]
  qfl <- qfl[pperiod]
  dataNew <- window(data, start=index(data)[warmup+1],end=index(data)[dim(data)[1]])
  dataNew$Q <- qfl
  data <- data.frame("datetime"=index(data)[pperiod],"P"=data[pperiod,"P"],"Q"=qfl)

  ## You can check out how this looks with the hydrograph function:
  # hydrograph(input=data.frame("datetime"=index(dataNew),"P"=coredata(dataNew$P),
  #                             "Q"=coredata(dataNew$Q)))

  ### STEP 2 ###
  # Identify discrete RR events from time series using hydromad::eventseq.
  # This function returns a zoo object, with core data consisting of an ordered
  # factor, representing the identified events, and the same time index as x.
  # Periods between events are left as NA, unless all = TRUE in which case they
  # are treated as separate events. The returned object stores thresh as an
  # attribute. Typical input parameters are:
  # x        = P and Q(quick) time series
  # thresh   = threshold value
  # inthresh = second threshold to define when events stop
  # mindur   = the minimum number of time steps in each event window
  # indur    = P must remain below inthresh for indur time steps in order to terminate an event
  # mingap   = the the minimum number of time steps that can separate events
  evpq <- eventseq(x       = dataNew[,c("P","Q")],
                  thresh   = c(0.1, median(qfl)),
                  inthresh = c(0,   median(qfl)),
                  mindur   = 6,
                  indur    = 1,
                  mingap   = 2)
  ## You can check out the structure of the evp object
  # str(evpq)
  #
  ## Plot
  xyplot(dataNew) + layer_(panel.xblocks(evpq, col = c("grey90", "grey80"), border = "grey80"))

  ### STEP 3 ###
  # Order events in descending order and calculate matching return periods
  eventsP <- sort(eventinfo(dataNew$P, evpq, FUN = max)$Value, decreasing = TRUE)
  returnPeriodP <- (length(eventsP)-1)/(1:length(eventsP))
  dP <- data.frame(x=returnPeriodP, y=eventsP)
  # loglogplot(dP)

  eventsQ <- sort(eventinfo(dataNew$Q, evpq, FUN = max)$Value, decreasing = TRUE)
  returnPeriodQ <- (length(eventsQ)-1)/(1:length(eventsQ))
  dQ <- data.frame(x=returnPeriodQ, y=eventsQ)
  # loglogplot(dQ)

  ### STEP 4 ###
  # Calculate the regression function to model dP and dQ

  # Calculate the regression function that describes dP
  # modelP.lm <- lm(formula = y ~ x + I(x^2), data = dP)
  modelP.lm <- lm(formula = y ~ log(x), data = dP)
  # Calculate the regression function that describes dQ
  # modelQ.lm <- lm(formula = y ~ x + I(x^2), data = dQ)
  modelQ.lm <- lm(formula = y ~ log(x), data = dQ)
  # Define some return periods
  Tr <- 2:round(max(returnPeriodP,returnPeriodQ),0)
  # Use predict to estimate the values for the return period.
  # Note that predict expects a data.frame and the col names need to match
  newP <- predict(modelP.lm, newdata = data.frame(x = Tr))
  newQ <- predict(modelQ.lm, newdata = data.frame(x = Tr))

  # plot(dP$x, dP$y, type="o", ylim=c(min(dP$y,newP),max(dP$y,newP)))
  # points(Tr, newP, col = "red")

  # plot(dQ$x, dQ$y, type="o", ylim=c(min(dQ$y,newQ),max(dQ$y,newQ)))
  # points(Tr, newQ, col = "red")

  # Final data.frame
  df <- data.frame("Tr"=Tr,"P"=newP,"Q"=newQ)

  return(df)

}
