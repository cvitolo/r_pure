PURE (R package)
================

Workflow for flood frequency analysis under uncertainty (work in progress!).

Project: PURE (Probability, Uncertainty and Risk in the Environment) - RACER consortium (Flood strand)

Check out: NERC-PURE programme at http://www.nerc.ac.uk/research/programmes/pure/.

Requirements:
sudo apt-get install libudunits2-dev
install.packages("hydroTSM")


#### Basics
Install and load packages
```R
# Install dependent packages from CRAN:
x <- c("zoo", "chron", "xts", "manipulate", "rgdal", 
       "sp", "gstat", "grid", "hydroTSM", "Hmisc", "raster", "reshape2", 
       "ggplot2", "qualV", "lhs", "MASS")
install.packages(x)
lapply(x, require, character.only=T); rm(x)

# Install dpendent package from R-Forge:
install.packages("fuse", repos="http://R-Forge.R-project.org")
library(fuse)

# Install dependent gists and packages from github:
library(devtools)
install_github("r_rnrfa", username = "cvitolo", subdir = "rnrfa")
library(rnrfa)

# Install pure package
install_github("r_pure", username = "cvitolo", subdir = "pure")
library(pure)
```

Make zoo objects for your time series: 

* Q is the streamflow time series. 
* E is the list of potential evapotranspiration time series. If E is unknown, wheather variables can be used to calculate the potential evapotranspiration (see section PET).
* P is the list of rainfall time series. If there are 3 raingauges in the catchment (e.g. P1, P2 and P3), the object P is: 
    + P <- list(P1,P2,P3)

An example is given below:
```R
# Load sample time series (this contains 3 objects: P, E and Q)
data(P1)
data(P2)
data(P3)
data(Q)
data(weather)
```

### Pre-processing
Below are a series of utility functions for time series pre-processing, divided in 4 categories: 

* Report
* Correct
* Aggregate
* Model specific preparation

##### Report
Report problems with time series using the function `ScanTS()`.
As an example you can use the example dataset provided with this package.

```R
# Report
ScanTS( P1, verbose = TRUE  )
ScanTS( P2, verbose = FALSE )
ScanTS( P3, verbose = FALSE )
ScanTS( Q,  verbose = FALSE )
ScanTS( weather$TD,  verbose = FALSE, returnNegInfo = FALSE )
ScanTS( weather$TW,  verbose = FALSE, returnNegInfo = FALSE )
ScanTS( weather$NR,  verbose = FALSE, returnNegInfo = FALSE )
ScanTS( weather$WS,  verbose = FALSE )
```

##### Correct, aggregate and prepare for modelling with FUSE
Often time series are recorded at non-regular time steps. You can shift your records to align them with a regular grid using the function `Irr2Reg()`.
```R
# From irregular to regular frequency time step:
P1Reg <- Irr2Reg( P1 )
P2Reg <- Irr2Reg( P2 )
P3Reg <- Irr2Reg( P3 )
QReg  <- Irr2Reg( Q )
TW    <- Irr2Reg( weather$TW )
TD    <- Irr2Reg( weather$TD )
NR    <- Irr2Reg( weather$NR )
WS    <- Irr2Reg( weather$WS )

# test the effect of Irr2Reg()
plot(P1[40:45])
lines(P1Reg[40:45],col="red")
```

Change any unrealistic values to NA (e.g. negative P and Q) using the function `CorrectNeg()`. 

```R
P1NoNeg <- CorrectNeg( P1Reg )

# test the effect of CorrectNeg()
plot(P1Reg)
lines(P1NoNeg,col="red")
```

Find coarser temporal resolution amongst a list of time series:
```R
myList <- list("P1" = P1NoNeg, "P2" = P2Reg, "P3" = P3Reg, 
               "Q" = QReg, "TW" = TW, "TD" = TD, "NR" = NR, "WS" = WS)
multiplier <- CommonTemporalResolution(myList); multiplier
```

Aggregate all the time series to the temporal resolution above:
```R
P1 <- aggregate(myList$P1, align.time(index(myList$P1), multiplier), FUN = sum)
P2 <- aggregate(myList$P2, align.time(index(myList$P2), multiplier), FUN = sum)
P3 <- aggregate(myList$P3, align.time(index(myList$P3), multiplier), FUN = sum)
TW <- aggregate(myList$TW, align.time(index(myList$TW), multiplier), FUN = mean)
TD <- aggregate(myList$TD, align.time(index(myList$TD), multiplier), FUN = mean)
NR <- aggregate(myList$NR, align.time(index(myList$NR), multiplier), FUN = mean)
WS <- aggregate(myList$WS, align.time(index(myList$WS), multiplier), FUN = mean)
Q  <- aggregate(myList$Q,  align.time(index(myList$Q ), multiplier), FUN = mean)

plot(myList$Q[1:100])
lines(Q[1:100],col="red")
```

Derive new variables, e.g. potential evapotranspiration from weather variables
```R
E <- pet(stationElevation=0,TD,TW,NR,WS)
```

Select periods with simultaneous recordings
```R
tsList <- list("P1" = P1, "P2" = P2, "P3" = P3, "E" = E, "Q" = Q)
newList <- ExtractOverlappingPeriod(tsList)
```

Aggregate in space, e.g. areal averaging using spatial interpolation methods
```R
tsList <- data.frame(index(newList),"P1"=newList$P1,
                     "P2"=newList$P2,"P2"=newList$P3)
P <- ArealAveraging(tsList,areas=c(0.3,0.6,0.1),interpolationMethod ="Thiessen")
```

Check if there are gaps in the records and infill
```R
any(is.na(P)) # FALSE

any(is.na(newList$E)) # This returns TRUE, therefore we will infill the missing values
Enomissing <- na.approx(newList$E)

any(is.na(newList$Q)) # This returns TRUE, therefore we will infill the missing values
Qnomissing <- na.approx(newList$Q)
```

If necessary, convert units to mm/day:
```R
P <- P*24 # from mm/h to mm/day
E <- Enomissing*24 # from mm/h to mm/day

Area <- 10.55 # Km2
Q <- Qnomissing*86.4/Area
```

Merge P, E and Q in 1 time series object
```R
DATA <- merge(P,E,Q)
```   

### Rainfall-Runoff modelling using FUSE
As an example, we could combine 50 parameter sets and 4 model structures to generate 200 model simulations.

Sample 50 parameter sets for FUSE, using LHS method
```R
library(fuse)
data(DATA)

set.seed(123)
NumberOfRuns <- 10    
parameters <- GeneratePsetsFUSE(NumberOfRuns)
```

Choose a list of models to take into account
```R
data(modlist)
parentModels <- c(60,230,342,426) # those are the parent models 
ModelList <- modlist[which(modlist$mid %in% parentModels),]
row.names(ModelList) <- NULL
```

Define the list of Model Performance Indices (MPIs)
```R
LAGTIME = function(x) lagtime(x$Qo,x$Qs)    
MAE     = function(x) mean(x$Qs - x$Qo, na.rm = TRUE)            
NSHF    = function(x) 1 - EF(x$Qo,x$Qs)           
NSLF    = function(x) 1 - EF( log(x$Qo) , log(x$Qs) )           
RR      = function(x) sum(x$Qs) /sum(x$Po)

MPIs <- list("LAGTIME"=LAGTIME,"MAE"=MAE,"NSHF"=NSHF,"NSLF"=NSLF,"RR"=RR)
```

Run simulations
```R
outputFolder <- "~"
deltim <- 1/24 # or multiplier/60/60/24
warmup <- round(dim(DATA)[1]/10,0)

MCsimulations(DATA,deltim,warmup,parameters,ModelList,outputFolder,MPIs)
```

### Find the best configuration(s) amongst those simulated
For this task it is necessary to install another experimental package called AMCA:

```R
# Install pure package
install_github("r_amca", username = "cvitolo", subdir = "amca")
library(amca)
```

Run the algorithm:
```R
results <- amca(DATA,ModelList,warmup,parameters,outputFolder)
```

The best configuration is stored in
```R
results$RETable
```

# Leave your feedback
I would greatly appreciate if you could leave your feedbacks via email (cvitolodev@gmail.com).
