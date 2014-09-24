PURE (R package)
================

Workflow for land management scenario analysis

Project: PURE (Probability, Uncertainty and Risk in the Environment) - RACER consortium (Flood strand)

Check out: NERC-PURE programme at http://www.nerc.ac.uk/research/programmes/pure/.

Requirements:
sudo apt-get install libudunits2-dev
install.packages("udunits2")
install.packages("outliers")
install.packages("hydroTSM")

### Basics
Install and load packages
```R
# Install dependent packages from CRAN:
x <- c("zoo", "chron", "xts", "manipulate", "udunits2", "outliers", "rgdal", 
       "sp", "gstat", "grid", "hydroTSM", "Hmisc", "raster", "reshape2", 
       "ggplot2", "lhs", "MASS")
install.packages(x)
lapply(x, require, character.only=T)
#library(plotly)
#source("~/Dropbox/Repos/github/r_uncertflood/makePlotly.R")

# Install dpendent package from R-Forge:
install.packages("fuse", repos="http://R-Forge.R-project.org")
library(fuse)

# Install dependent gists and packages from github:
library(devtools)

install_github("r_amca", username = "cvitolo", subdir = "amca")
library(amca)

install_github("r_rnrfa", username = "cvitolo", subdir = "rnrfa")
library(rnrfa)
source_gist("https://gist.github.com/cvitolo/f9d12402956b88935c38")

# Install pure package
install_github("r_pure", username = "cvitolo", subdir = "pure")
library(pure)
```

### Use PURE data
Use data from Imperial College database:
```R
CatchmentName <- "Pontbren" 
SubcatchmentName <- 9      
deltim <- 1/24
# do not use ~ for home folder, rgdal does not like it!
datafolder <- "/home/claudia/Dropbox/Projects/PURE/PURE_shared/Data/"
```
Load Time series data and GIS layers
```R
DataList <- LoadMyData(CatchmentName,SubcatchmentName,datafolder)
# manipulatePlot(DataList)
```

# Pre-processing: Utility functions for Time Series screening

### Scan the time series and report problems
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

### Convert from irregular to regular time series
This function shifts the records to align with a regular grid
```R
regTS <- irreg2regTS(CatchmentName, DataList, deltim)
```

### Check if there are gaps in the records and infill 
```R
gaps <- findGaps(regTS,deltim)
NoGaps <- fullrangeTS(regTS, gaps$fullranges)
infilled <- fillGaps(NoGaps)
``` 

### Correct negative values (if applicable)
```R
NoNeg <- correctNeg(infilled)
```

### ToDo: Remove Outliers?

### Areal averaging
Available algorithms: "AritmeticMean", "Thiessen", "IDW", "OK" 
```R
InterpolationMethod <- "Thiessen"
Paveraged <- RainfallArealAveraging(CatchmentName,overlappingTS(NoNeg)$P,DataList,InterpolationMethod) 
```

### Convertion to use TS with FUSE
```R
Pconverted <- Paveraged # original values are already in mm/d
Econverted <- ConvertTS(NoNeg$E,from="mm/h", to="mm/d")
Qconverted <- ConvertTS(NoNeg$Q,from="l/s", to="mm/day",
                        optionalInput=DataList$Area*10^12)
```

### Extract overlapping periods
```R
DATA <- ExtractOverlappingPeriod(Pconverted,Econverted,Qconverted,ignoreQ=TRUE)
```

# Rainfall-Runoff modelling with FUSE

### Generate parameter set to use in FUSE
```R
set.seed(123)
NumberOfRuns <- 1000
parameters <- GeneratePsetsFUSE(NumberOfRuns)
```

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
