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


#### Basics
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

### Pre-processing
Below are a series of utility functions for time series screening. 

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

##### Convert from irregular to regular time series
This function shifts the records to align with a regular grid
```R
regTS <- irreg2regTS(CatchmentName, DataList, deltim)
```

##### Check if there are gaps in the records and infill 
```R
gaps <- findGaps(regTS,deltim)
NoGaps <- fullrangeTS(regTS, gaps$fullranges)
infilled <- fillGaps(NoGaps)
``` 

##### Correct negative values (if applicable)
```R
NoNeg <- correctNeg(infilled)
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
NumberOfRuns <- 1000
parameters <- GeneratePsetsFUSE(NumberOfRuns)
```

