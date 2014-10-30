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
plot(DATA)

### Rainfall-Runoff modelling using FUSE
# As an example, we could combine 50 parameter sets and 4 model structures to generate 200 model simulations.

# Sample 50 parameter sets for FUSE, using LHS method
library(fuse)

set.seed(123)
NumberOfRuns <- 50
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
results <- amca(DATA,
                parameters,
                MPIs,
                outputFolder,
                selectedModels = ModelList$mid,
                warmup,
                verbose=TRUE)

# The best configuration is stored in results$RETable

# Visualise ensemble results
PlotEnsembles(results$BoundsIE$bounds, results$BoundsRE$discharges,
              label1="IE min-max", label2="RE percentiles")

PlotModelSimilarities(ModelList,results$RETable)
#PlotParameterSimilarities(results$RETable, parameters)

save(DATA,parameters,MPIs,ModelList,warmup, results, file="~/testPontbren.rda")

###########################################
# Frequency Analysis of  TIME SERIES DATA #
###########################################
load("~/testPontbren.rda")
library(pure)
library(hydromad)

# return statistics to identify threshold
source('~/Dropbox/Repos/github/r_pure/pure/R/EventIdentification.R')
plot(DATA)
summary(DATA$P)
summary(DATA$Q)

df <- EventIdentification(DATA)

################
# Curve Number #
################
source('~/Dropbox/Repos/github/r_pure/pure/R/calculateCN.R')

# 3. Determine the CN for each event
# where P & Q are in inches and area is in acre
Q <- df$Q/25.4
P <- df$P/25.5
area <- 4.8*247.105 # area <- DataList$Area*247.105

CN <- calculateCN(P,Q); df$CN <- CN
myCN <- median( sort(CN, decreasing = TRUE)[1:5] )

