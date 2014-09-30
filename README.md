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

# install_github("r_amca", username = "cvitolo", subdir = "amca")
library(amca)

# install_github("r_rnrfa", username = "cvitolo", subdir = "rnrfa")
library(rnrfa)
source_gist("https://gist.github.com/cvitolo/f9d12402956b88935c38")

# Install pure package
# install_github("r_pure", username = "cvitolo", subdir = "pure")
library(pure)

# library(plotly)
# source("~/Dropbox/Repos/github/r_uncertflood/makePlotly.R")
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
```

### Pre-processing
Below are a series of utility functions for time series pre-processing. Those are divided into 4 sections: 

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
ScanTS( E,  verbose = FALSE, returnNegInfo = FALSE )
ScanTS( Q,  verbose = FALSE )
```

##### Correct
Correct unrealistis values (e.g. negative P and Q) using the function `CorrectNeg()`. 
Note E can be negative!
```R
P1NoNeg <- CorrectNeg( P1 )
```

Often time series are recorded at non-regular time steps. You can shift your records to align them with a regular grid using the function `Irr2Reg()`.
```R
# set a regular time step in days:
P1Reg <- Irr2Reg( P1NoNeg )
P2Reg <- Irr2Reg( P2 )
P3Reg <- Irr2Reg( P3 )
EReg  <- Irr2Reg( E )
QReg  <- Irr2Reg( Q )

plot(P1NoNeg[40:45])
lines(P1Reg[40:45],col="red")
```

# Leave your feedback
I would greatly appreciate if you could leave your feedbacks via email (cvitolodev@gmail.com).
