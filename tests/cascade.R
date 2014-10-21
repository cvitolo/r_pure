# Install/load libraries
x <- c("zoo", "chron", "xts", "manipulate", "rgdal", "tgp", "rgdal",
       "sp", "gstat", "grid", "hydroTSM", "Hmisc", "raster", "reshape2",
       "ggplot2", "qualV", "lhs", "MASS",
       "amca","fuse","pure")
# install.packages(x)
lapply(x, require, character.only=T); rm(x)
source('~/Dropbox/Repos/github/r_pure/tests/LoadPureData.R')

# set the path to the input files folder (NOTE: do not use ~ for home folder,
# rgdal does not like it!)
datafolder <- "/home/claudia/Dropbox/Projects/PURE/PURE_shared/Data/"

# Extract data for Pontbren, subcatchment 9
CatchmentName <- "Pontbren"
SubcatchmentName <- 9
deltim <- 1/24
DataList <- LoadPureData(CatchmentName,SubcatchmentName,datafolder)

# explore the content of DataList using:
# names(DataList)
# str(DataList)
# in this case there are 3 rain gauges in the area: "P2" "P3" "P5" (10 minutes)
# evapotranspiration calculated from MORECS (1 time series, 1 hour)
# streamflow (1 time series, 15 minutes)

# Report
ScanTS( DataList$P$P2, verbose = TRUE  )
ScanTS( DataList$P$P3)
ScanTS( DataList$P$P5)
ScanTS( DataList$E, returnNegInfo = FALSE )
ScanTS( DataList$Q )

##### Correct, aggregate and prepare for modelling with FUSE
# Often time series are recorded at non-regular time steps.
# You can shift your records to align them with a regular grid using the
# function `Irr2Reg()`.

# From irregular to regular frequency time step:
P2Reg <- Irr2Reg( DataList$P$P2 )
P3Reg <- Irr2Reg( DataList$P$P3 )
P5Reg <- Irr2Reg( DataList$P$P5 )
QReg  <- Irr2Reg( DataList$Q )

# test the effect of Irr2Reg()
plot(DataList$P$P2[1:100])
lines(P2Reg[1:100],col="red")

# Pontbren does not have any unrealistic values.
# If there were some, they could have been changed to NA
# using the function `CorrectNeg()`.

# Find coarser temporal resolution (in seconds!) amongst a list of time series
x <- CommonTemporalResolution(list("P2" = P2Reg,
                                   "P3" = P3Reg,
                                   "P5" = P5Reg,
                                   "E"  = DataList$E,
                                   "Q"  = QReg)      )

# and aggregate all of them to the same temporal resolution: 1 hour
# P and E should be aggregated with FUN = mean for rates (FUN = sum for volumes)
# Q should always be aggregated with FUN = mean
P2a <- period.apply(P2Reg,      endpoints(P2Reg,      "seconds", x), FUN = mean)
P3a <- period.apply(P3Reg,      endpoints(P3Reg,      "seconds", x), FUN = mean)
P5a <- period.apply(P5Reg,      endpoints(P5Reg,      "seconds", x), FUN = mean)
Ea  <- period.apply(DataList$E, endpoints(DataList$E, "seconds", x), FUN = mean)
Qa  <- period.apply(QReg,       endpoints(QReg,       "seconds", x), FUN = mean)

# format time index consistently
P2a_align <- align.time(P2a, x)
P3a_align <- align.time(P3a, x)
P5a_align <- align.time(P5a, x)
Ea_align  <- align.time(xts(Ea),  x)
Qa_align  <- align.time(Qa,  x)

# test the effect of the aggregation
plot(P2Reg[1:100])
lines(P2a_align[1:100], col="red")

# There are no new variables to derive, e.g. potential evapotranspiration
# from weather variables (in case, use the function pet())

# Select periods with simultaneous recordings
newList <- CommonRecordingTime(list("P2" = P2a_align,
                                    "P3" = P3a_align,
                                    "P5" = P5a_align,
                                    "E"  = Ea_align,
                                    "Q"  = Qa_align) )

# Aggregate in space, e.g. areal averaging using spatial interpolation methods
tsList <- data.frame(index(newList$P2),"P2"=newList$P2,
                     "P3"=newList$P3,"P5"=newList$P5)
P <- ArealAveraging(tsList,  areas=DataList$A, interpolationMethod ="Thiessen")

# Check if there are gaps in the records and infill
any(is.na(P))         # FALSE
any(is.na(newList$E)) # FALSE
any(is.na(newList$Q)) # FALSE

# If necessary, convert units to mm/day:
# P was already recorded in mm/d, no need for convertion
E <- newList$E*24          # from mm/h to mm/day
Q <- newList$Q*0.0864/DataList$Area # from l/s to mm/day (area is in Km2)

# Merge P, E and Q in 1 time series object
DATA <- merge(P,E,Q); names(DATA) <- c("P","E","Q")

### Rainfall-Runoff modelling using FUSE
# As an example, we could combine 50 parameter sets and 4 model structures to generate 200 model simulations.

# Sample 50 parameter sets for FUSE, using LHS method
library(fuse)
data(DATA)

set.seed(123)
NumberOfRuns <- 10
parameters <- GeneratePsetsFUSE(NumberOfRuns)

# Choose a list of models to take into account
data(modlist)
parentModels <- c(60,230,342,426) # those are the parent models
ModelList <- modlist[which(modlist$mid %in% parentModels),]
row.names(ModelList) <- NULL

# Define the list of Model Performance Indices (MPIs)
library(tiger)
library(qualV)

LAGTIME = function(x) lagtime(x$Qo,x$Qs)
MAE     = function(x) mean(x$Qs - x$Qo, na.rm = TRUE)
NSHF    = function(x) 1 - EF(x$Qo,x$Qs)
NSLF    = function(x) 1 - EF( log(x$Qo) , log(x$Qs) )
RR      = function(x) sum(x$Qs) /sum(x$Po)

MPIs <- list("LAGTIME"=LAGTIME,"MAE"=MAE,"NSHF"=NSHF,"NSLF"=NSLF,"RR"=RR)

# Run simulations
outputFolder <- "~"
deltim <- 1/24 # or dt/60/60/24
warmup <- round(dim(DATA)[1]/10,0)

# It is recommended to run simulations on HPC facilities.
# However small batches can be run locally using the function MCsimulations()
MCsimulations(DATA,deltim,warmup,parameters,ModelList,outputFolder,MPIs)


### Find the best configuration(s) amongst those simulated

# Run the algorithm
results <- amca(DATA,ModelList,warmup,parameters,outputFolder)

# The best configuration is stored in
results$RETable

























exampleTimeSeries <- list("P"=list("P1"=DataList$P$P2,"P2"=DataList$P$P3,"P3"=DataList$P$P5),"E"=DataList$E,"Q"=DataList$Q)
save(exampleTimeSeries,file="data/exampleTimeseries.rda")

# manipulatePlot(DataList)

##################
# Pre-processing #
##################

# ########## Utility functions for Time Series ##########

#ScanTS( DataList$P, returnGapsInfo=TRUE, returnTimeInfo=TRUE, returnNegInfo=TRUE )  # units = mm/d
#ScanTS( DataList$E, returnGapsInfo=TRUE, returnTimeInfo=TRUE, returnNegInfo=FALSE ) # units = mm/h
#ScanTS( DataList$Q, returnGapsInfo=TRUE, returnTimeInfo=TRUE, returnNegInfo=TRUE )  # units = l/s

# ########## Pre-processing ##########

# Convertion from irregular to regular time series (this function shifts the records to align with a regular grid)
regTS <- irreg2regTS(CatchmentName, DataList, deltim)
# manipulatePlot(regTS)
# the two plots look different, we can check here:
# plot(regTS$P[[1]],type="l")
# lines(DataList$P[[1]],col="red",lty=2)
# The new time series reach higher values, this is normal as we cumulated the value of subtimesteps.
# if we check the difference in volumes, we find that there is no error.
# sum(DataList$P[[1]])-sum(regTS$P[[1]])

# Check if there are gaps in the records and infill
gaps <- findGaps(regTS,deltim)
NoGaps <- fullrangeTS(regTS,gaps$fullranges)
# manipulatePlot(NoGaps)

infilled <- fillGaps(NoGaps)
# manipulatePlot(infilled)

# Negative values to be corrected (if applicable)
NoNeg <- correctNeg(infilled)

# ToDo: Remove Outliers?

# RAINFALL AREAL AVERAGING
# Available algorithms: "AritmeticMean", "Thiessen", "IDW", "OK"
InterpolationMethod <- "Thiessen"
Paveraged <- RainfallArealAveraging(CatchmentName,overlappingTS(NoNeg)$P,DataList,InterpolationMethod)

# Convertion to use TS with FUSE
# Pconverted <- ConvertTS(Paveraged,from="mm/h", to="mm/d") # original values are already in mm/d
Pconverted <- Paveraged
Econverted <- ConvertTS(NoNeg$E,from="mm/h", to="mm/d")
Qconverted <- ConvertTS(NoNeg$Q,from="l/s", to="mm/day",optionalInput=DataList$Area*10^12)

# Extract overlapping periods
DATA <- ExtractOverlappingPeriod(Pconverted,Econverted,Qconverted,ignoreQ=TRUE)

#######################################
# Rainfall-Runoff modelling with FUSE #
#######################################

# Generate parameter set to use in FUSE
set.seed(123)
NumberOfRuns <- 50    #1000
parameters <- GeneratePsetsFUSE("NumberOfRuns"=NumberOfRuns,SamplingType="LHS",RoutingParameterRange=c(0.01,5))

#parentModels <- c(60,230,342,426) # 60=topmodel, 230 = arno-vic, 342 = prms, 426 = sacramento

# 4 figures arranged in 2 rows and 2 columns
topmodel <- FUSEcalibration("DATA"=DATA, "deltim"=deltim, "mid"=60, "parameters"=parameters)
arnovic <- FUSEcalibration("DATA"=DATA, "deltim"=deltim, "mid"=230, "parameters"=parameters)
prms <- FUSEcalibration("DATA"=DATA, "deltim"=deltim, "mid"=342, "parameters"=parameters)
sacramento <- FUSEcalibration("DATA"=DATA, "deltim"=deltim, "mid"=426, "parameters"=parameters)

# multiplot(topmodel$p, arnovic$p, prms$p, sacramento$p, cols=2)

# compareModels <- function(column=2){
#   plot(topmodel$t[,column],type="l")
#   lines(arnovic$t[,column],col="red")
#   lines(prms$t[,column],col="blue")
#   lines(sacramento$t[,column],col="green")
# }
# manipulate(compareModels(column),column=slider(2,4))

#save(topmodel,arnovic,prms,sacramento,file="~/Dropbox/Repos/r_uncertflood/parentmodelsimulations.rda")



################
# Curve Number #
################
library(hydromad)
### Frequency Analysis of  TIME SERIES DATA

# 1. Partition the rainfall-runoff time series into individual events and calculate rainfall and storm runoff depths [Boorman et al., 1995]
x <- DATA[,c("P","Q")]
evp <- eventseq(x$P, thresh = 5, inthresh = 1, indur = 4, continue = TRUE)
evq <- eventseq(x$Q, thresh = 2, indur = 4, mingap = 5)

# 2. Order the data, i.e., sort the individual rainfall and runoff depths independently in descending order to match rainfall and runoff return periods
eventsQ <- sort(eventinfo(x$Q, evq)$Value, decreasing = TRUE)
returnPeriodQ <- (length(eventsQ)-1)/(1:length(eventsQ))
# exceedenceProbability <- 1/returnPeriod
dQ <- data.frame(x=returnPeriodQ, y=eventsQ)
source('~/Dropbox/Repos/github/r_uncertflood/loglogplot.R')
loglogplot(dQ)

eventsP <- sort(eventinfo(x$P, evp)$Value, decreasing = TRUE)
returnPeriodP <- (length(eventsP)-1)/(1:length(eventsP))
# exceedenceProbability <- 1/returnPeriod
dP <- data.frame(x=returnPeriodP, y=eventsP)
loglogplot(dP)

model.lm <- lm(formula = y ~ x + I(x^2) + I(x^3) + I(x^4), data = dP)
# Use predict to estimate the values for the return period.
# Note that predict expects a data.frame and the col names need to match
newY <- predict(model.lm, newdata = data.frame(x = dQ$x))

plot(dP$x, dP$y, type="o")
points(dQ$x, newY, col = "red")

# Final data.frame
df <- data.frame(Tr=dQ$x,P=newY,Q=dQ$y)

# 3. Determine the CN for each event
# where P & Q are in inches and area is in acre
Q <- df$Q/25.4
P <- df$P/25.5
area <- DataList$Area*247.105

source('~/Dropbox/Repos/github/r_uncertflood/UncertFlood/R/calculateCN.R')
CN <- calculateCN(P,Q)
df$CN <- CN

result <- nls(log(CN)~log(a + (100 - a) * exp(-b*P)),start=list(a=1,b=1),data=df)
summary(result)

a <- summary(result)$coefficients["a","Estimate"] # CN*USDA
b <- summary(result)$coefficients["b","Estimate"] # k

### SPATIAL DATA

catchment <- DataList$catchment
DTM <- DataList$DTM
Soil <- DataList$Soil
LandUse <- DataList$LandUse

plot(Soil)
plot(catchment,col=rgb(1,0,0,0.2),add=TRUE)

# Mapping HOST to USDA classes
t <- data.frame(table(raster::extract(Soil,catchment)))
t <- cbind(t,"USDAclass"=rep(NA, dim(t)[1]))
for (r in 1: dim(t)[1]) {
  if (any(c(1,2,3,5,11,13)==t[r,1])) t[r,3] <- "A"
  if (any(c(4,7)==t[r,1])) t[r,3] <- "AB"
  if (any(c(6,8,9,10,16)==t[r,1])) t[r,3] <- "B"
  if (17==t[r,1]) t[r,3] <- "BC"
  if (any(c(18,19,20)==t[r,1])) t[r,3] <- "C"
  if (any(c(14,15,28)==t[r,1])) t[r,3] <- "CD"
  if (any(c(12,21,22,23,24,25,26,27,29)==t[r,1])) t[r,3] <- "D"
}

USDAclass <- unique(t[,3])
taggregated <- data.frame(USDAclass,"Ncells"=rep(NA,length(USDAclass)))
for (r in 1:length(USDAclass)) {
  taggregated[r,2] <- sum(t[which(t[,3]==taggregated[r,1]),2])
}

plot(LandUse)
plot(catchment,col=rgb(1,0,0,0.2),add=TRUE)

# Mapping Land Cover Map to USDA classes
t2 <- data.frame(table(extract(LandUse,catchment)))


################################################################################

##### Scan the time series and report problems
```R
# Precipitation units = mm/d
ScanTS(DataList$P, returnGapsInfo=TRUE,
       returnTimeInfo=TRUE, returnNegInfo=TRUE)

# Potential evapotranspiration units = mm/h
ScanTS(DataList$E, returnGapsInfo=TRUE,
       returnTimeInfo=TRUE, returnNegInfo=FALSE)

# Streamflow discharge units = l/s
ScanTS(DataList$Q, returnGapsInfo=TRUE,
       returnTimeInfo=TRUE, returnNegInfo=TRUE)
```

##### Correct unrealistis values (e.g. negative P and Q)
```R
NoNeg <- correctNeg(infilled)
```

##### Correct unrealistis values (e.g. negative P and Q)
```R
NoNeg <- correctNeg(infilled)
```

##### Convert from irregular to regular time series
This function shifts the records to align them with a regular grid
```R
regTS <- irreg2regTS(CatchmentName, DataList, deltim)
```

##### Check if there are gaps in the records and infill
```R
gaps <- findGaps(regTS,deltim)
NoGaps <- fullrangeTS(regTS, gaps$fullranges)
infilled <- fillGaps(NoGaps)
```



##### ToDo: Remove Outliers?

##### Areal averaging
Available algorithms: "AritmeticMean", "Thiessen", "IDW", "OK"
```R
InterpolationMethod <- "Thiessen"
Paveraged <- RainfallArealAveraging(CatchmentName,overlappingTS(NoNeg)$P,DataList,InterpolationMethod)
```

##### Convertion to use TS with FUSE
```R
Pconverted <- Paveraged # original values are already in mm/d
Econverted <- ConvertTS(NoNeg$E,from="mm/h", to="mm/d")
Qconverted <- ConvertTS(NoNeg$Q,from="l/s", to="mm/day",
                        optionalInput=DataList$Area*10^12)
```

##### Extract overlapping periods
```R
DATA <- ExtractOverlappingPeriod(Pconverted,Econverted,Qconverted,ignoreQ=TRUE)
```

# Rainfall-Runoff modelling with FUSE

##### Generate parameter set to use in FUSE
```R
set.seed(123)
NumberOfRuns <- 50
parameters <- GeneratePsetsFUSE(NumberOfRuns)
```

##### Ensemble example usage
Define a group of model structures to use
```R
mids <- c(60, 230, 342, 426)
```

Run a multi-model calibration using the Nash-Sutcliffe efficiency as objective function
```R
indices <- rep(NA,4*NumberOfRuns)
discharges <- matrix(NA,ncol=4*NumberOfRuns,nrow=dim(DATA)[1])
kCounter <- 0

for (m in 1:4){

  myMID <- mids[m]

  for (pid in 1:NumberOfRuns){

    kCounter <- kCounter + 1
    ParameterSet <- as.list(parameters[pid,])

    Qrout <- RunFUSE(DATA, parameters[pid,], deltim, myMID)

    indices[kCounter] <- EF(DATA$Q,Qrout)
    discharges[,kCounter] <- Qrout

  }
}
```

Compare results
```R
bestRun <- which(indices == max(indices))

bestModel <- function(runNumber){
  if (runNumber<(NumberOfRuns+1)) myBestModel <- "TOPMODEL"
  if (runNumber>(NumberOfRuns+1) & runNumber<(2*NumberOfRuns+1)) myBestModel <- "ARNOXVIC"
  if (runNumber>(2*NumberOfRuns+1) & runNumber<(3*NumberOfRuns+1)) myBestModel <- "PRMS"
  if (runNumber>(3*NumberOfRuns+1) & runNumber<(4*NumberOfRuns+1)) myBestModel <- "SACRAMENTO"
  return(myBestModel)
}
bestModel(bestRun)

plot(coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=0.5)

for(pid in 1:(4*NumberOfRuns)){
  lines(discharges[,pid], col="gray", lwd=0.1)
}

lines(coredata(DATA$Q),col="black", lwd=1)
lines(discharges[,bestRun],col="red", lwd=1)
```

How the best simulations of each model structure compare to each other?
```R
bestRun0060 <- which(indices[1:NumberOfRuns] == max(indices[1:NumberOfRuns]))
bestRun0230 <- NumberOfRuns + which(indices[(NumberOfRuns+1):(2*NumberOfRuns)] == max(indices[(NumberOfRuns+1):(2*NumberOfRuns)]))
bestRun0342 <- 2*NumberOfRuns + which(indices[(2*NumberOfRuns+1):(3*NumberOfRuns)] == max(indices[(2*NumberOfRuns+1):(3*NumberOfRuns)]))
bestRun0426 <- 3*NumberOfRuns + which(indices[(3*NumberOfRuns+1):(4*NumberOfRuns)] == max(indices[(3*NumberOfRuns+1):(4*NumberOfRuns)]))

plot(coredata(DATA$Q),type="l",xlab="",ylab="Streamflow [mm/day]", lwd=1)
lines(discharges[,bestRun0060], col="green", lwd=1)
lines(discharges[,bestRun0230], col="blue", lwd=1)
lines(discharges[,bestRun0342], col="pink", lwd=1)
lines(discharges[,bestRun0426], col="orange", lwd=1)

legend("top", horiz=TRUE, cex=0.65,
       c("TOPMODEL", "ARNOXVIC", "PRMS","SACRAMENTO"),
       col = c("green", "blue", "pink", "orange"),
       lty = c(1, 1, 1, 1))
```

##### GLUE UNCERTAINTY ANALYSIS
glue()
