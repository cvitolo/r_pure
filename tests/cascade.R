

# install_github("r_amca", username = "cvitolo", subdir = "amca")
# library(amca)



# library(plotly)
# source("~/Dropbox/Repos/github/r_uncertflood/makePlotly.R")

#for Nathan's functions
library(lhs)
library(MASS)

#for Claudia's functions
library(manipulate)
library(Hmisc)
library(ggplot2)
library(reshape)
library(raster)
library(xts)
#library(plotly)
#source("~/Dropbox/Repos/github/r_uncertflood/makePlotly.R")

library(fuse)
library(damach)
library(uncertflood)
library(rnrfa)
# set bounding box
bb <- list("LonMin"=-3.4477329, "LonMax"=-3.3914280, "LatMin"=52.632490, "LatMax"=52.656863)
x <- getStationSummary(bb$LonMin, bb$LonMax, bb$LatMin, bb$LatMax)

# Use data from PURE database:
CatchmentName <- "Pontbren"
SubcatchmentName <- 9
deltim <- 1/24
# do not use ~ for home folder, rgdal does not like it!
datafolder <- "/home/claudia/Dropbox/Projects/PURE/PURE_shared/Data/"
# datafolder <- ifelse(getwd()=="/home/cvitolo/ModelCascade","/home/cvitolo/ModelCascade/data/","/home/claudia/Dropbox/PURE/Data/")

# ########## Load data (TS + GIS files) ##########

DataList <- LoadMyData(CatchmentName,SubcatchmentName,datafolder)

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
