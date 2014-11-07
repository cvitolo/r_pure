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
x <- c("zoo", "chron", "xts", "manipulate", "rgdal", "tgp",
       "sp", "gstat", "grid", "hydroTSM", "Hmisc", "raster", "reshape2", 
       "ggplot2", "qualV", "lhs", "MASS")
install.packages(x)
lapply(x, require, character.only=T); rm(x)

# Install dpendent package from R-Forge:
install.packages("fuse", repos="http://R-Forge.R-project.org")

# Install dependent gists and packages from github:
library(devtools)
install_github("cvitolo/r_rnrfa", subdir = "rnrfa")

# Install pure package
install_github("cvitolo/r_pure", subdir = "pure")
```

Make zoo objects for your time series: 

* Q is the streamflow time series. 
* E is the list of potential evapotranspiration time series. If E is unknown, wheather variables can be used to calculate the potential evapotranspiration (see section PET).
* P is the list of rainfall time series. If there are 3 raingauges in the catchment (e.g. P1, P2 and P3), the object P is: 
    + P <- list(P1,P2,P3)

An example is given below:
```R
library(pure)
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
```

##### Common temporal resolution
Find coarser temporal resolution (in seconds!) amongst a list of time series and aggregate all of them to the same temporal resolution
```R
myList <- list("P1" = P1NoNeg, "P2" = P2Reg, "P3" = P3Reg, 
               "Q" = QReg, "TW" = TW, "TD" = TD, "NR" = NR, "WS" = WS)
x <- CommonTemporalResolution(myList)

# use "mean" for rates, "sum" for volumes
P1 <- period.apply(myList$P1, endpoints(myList$P1, "seconds", x), mean)
P2 <- period.apply(myList$P2, endpoints(myList$P2, "seconds", x), mean)
P3 <- period.apply(myList$P3, endpoints(myList$P3, "seconds", x), mean)
TW <- period.apply(myList$TW, endpoints(myList$TW, "seconds", x), mean)
TD <- period.apply(myList$TD, endpoints(myList$TD, "seconds", x), mean)
NR <- period.apply(myList$NR, endpoints(myList$NR, "seconds", x), mean)
WS <- period.apply(myList$WS, endpoints(myList$WS, "seconds", x), mean)
Q  <- period.apply(myList$Q,  endpoints(myList$Q,  "seconds", x), mean)
```

Derive new variables, e.g. potential evapotranspiration from weather variables
```R
E <- pet(stationElevation=0,TD,TW,NR,WS)
```

Align time stamps
```R
# format time index consistently
P1a <- align.time(P1, x)
P2a <- align.time(P2, x)
P3a <- align.time(P3, x)
Ea <- align.time(E, x)
Qa <- align.time(Q, x)

```

Select periods with simultaneous recordings
```R
# Select periods with simultaneous recordings
newList <- CommonRecordingTime(list("P1" = P1a, "P2" = P2a, "P3" = P3a, 
                                    "E" = Ea, "Q" = Qa) )
```

Aggregate in space, e.g. areal averaging using spatial interpolation methods
```R
tsList <- data.frame(index(newList),"P1"=newList$P1,
                     "P2"=newList$P2,"P2"=newList$P3)
                     
P <- ArealAveraging(tsList,areas=c(0.3,0.6,0.1),
                    interpolationMethod ="Thiessen")
```

Check if there are gaps in the records and infill
```R
any(is.na(P)) # FALSE

any(is.na(newList$E)) # This returns TRUE, infill the missing values
Enomissing <- na.approx(newList$E)

any(is.na(newList$Q)) # This returns TRUE, infill the missing values
Qnomissing <- na.approx(newList$Q)
```

If necessary, convert units to mm/day:
```R
# P is already in mm/day
E <- Enomissing*24 # from mm/h to mm/day

Area <- 10.55 # Km2
Q <- Qnomissing*0.086.4/Area # from l/s to mm/day
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
parentModels <- c(60,230,342,426) # those are the parent models
```

Define the list of Model Performance Indices (MPIs)
```R
library(tiger)
library(qualV)

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
deltim <- 1/24 # or dt/60/60/24
warmup <- round(dim(DATA)[1]/10,0)

# It is recommended to run simulations on HPC facilities. 
# However small batches can be run locally using the function MCsimulations()
MCsimulations(DATA,deltim,warmup,parameters,parentModels,outputFolder,MPIs)
```

### Find the best configuration(s) amongst those simulated
For this task it is necessary to install another experimental package called AMCA:

```R
# Install pure package
install_github("cvitolo/r_amca", subdir = "amca")
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

### Frequency Analysis
For this task it is necessary to install another packages:

```R
# Install and load packages 
install.packages("hydromad", repos="http://hydromad.catchment.org")
install.packages("EcoHydRology")
library(EcoHydRology)
```

Identify events Precipitation-discharge events:
```R
df <- EventIdentification(DATA)
```

### Curve Number
Determine the CN for each event where P & Q are in inches and area is in acre
```R
Q <- df$Q/25.4       # in inches
P <- df$P/25.5       # in inches
area <- Area*247.105 # in acre

CN <- calculateCN(P,Q); df$CN <- CN
myCN <- median( sort(CN, decreasing = TRUE)[1:5] )

```

# Leave your feedback
I would greatly appreciate if you could leave your feedbacks via email (cvitolodev@gmail.com).
