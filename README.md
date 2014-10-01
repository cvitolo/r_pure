PURE (R package)
================

Workflow for flood frequency analysis under uncertainty.

Project: PURE (Probability, Uncertainty and Risk in the Environment) - RACER consortium (Flood strand)

Check out: NERC-PURE programme at http://www.nerc.ac.uk/research/programmes/pure/.

Requirements:
sudo apt-get install libudunits2-dev
install.packages("udunits2")
install.packages("outliers")
install.packages("hydroTSM")


#### Basics
Install and load packages
```R
# Install dependent packages from CRAN:
x <- c("zoo", "chron", "xts", "manipulate", "udunits2", "outliers", "rgdal", 
       "sp", "gstat", "grid", "hydroTSM", "Hmisc", "raster", "reshape2", 
       "ggplot2", "qualV", "lhs", "MASS")
# install.packages(x)
lapply(x, require, character.only=T)

# Install dpendent package from R-Forge:
# install.packages("fuse", repos="http://R-Forge.R-project.org")
library(fuse)

# Install dependent gists and packages from github:
library(devtools)

# install_github("r_rnrfa", username = "cvitolo", subdir = "rnrfa")
library(rnrfa)

# Install pure package
# install_github("r_pure", username = "cvitolo", subdir = "pure")
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
data(E)
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

plot(P1[40:45])
lines(P1Reg[40:45],col="red")
```

Change any unrealistic values to NA (e.g. negative P and Q) using the function `CorrectNeg()`. 

```R
P1NoNeg <- CorrectNeg( P1Reg )

plot(P1Reg)
lines(P1NoNeg,col="red")
```

Set a common temporal resolution and aggregate all the time series
```R
myList <- list("P1" = P1NoNeg, "P2" = P2Reg, "P3" = P3Reg, 
               "E" = EReg, "Q" = QReg, 
               "TW" = TW, "TD" = TD, "NR" = NR, "WS" = WS)
multiplier <- CommonTemporalResolution(myList)

temp <- index(as.xts(aggregate(P1, align.time(index(P1), multiplier))))
P <- zoo(aggregate(P1NoNeg, mean))

P <- aggregate(P1NoNeg,by=rep(temp,), FUN=sum)
```

Derive new variables, e.g. E from weather variables
```R
E <- pet(stationElevation=0,TD,TW,NR,WS)
```

Aggregate in space, e.g. averal averaging using spatial interpolation methods

Check if there are gaps in the records and infill
```R
gaps <- findGaps(regTS,deltim)
NoGaps <- fullrangeTS(regTS, gaps$fullranges)
infilled <- fillGaps(NoGaps)
```

Convert units

##### Aggregate and convert based on FUSE model requirements


# Leave your feedback
I would greatly appreciate if you could leave your feedbacks via email (cvitolodev@gmail.com).
