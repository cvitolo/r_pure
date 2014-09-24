#' Scans time series given in input and checks for un-even recording times,
#' negative values and missing values.
#'
#' @param OriginalDataset Time series, it can be zoo, xts, etc
#' @param returnGapsInfo bolean, default value is TRUE
#' @param returnTimeInfo bolean, default value is TRUE
#' @param returnNegInfo bolean, default value is TRUE
#' @param verbose bolean, default value is TRUE
#'
#' @return print useful information
#'
#' @examples
#' # ScanTS( dataset=Q[[8]] , returnGapsInfo=TRUE, returnTimeInfo=TRUE,
#' #         returnNegInfo=TRUE , verbose = TRUE)
#'

ScanTS <- function(OriginalDataset,
                   returnGapsInfo=TRUE,
                   returnTimeInfo=TRUE,
                   returnNegInfo=TRUE,
                   verbose = TRUE){

  ListOfDatasets <- as.list(OriginalDataset)

  for (counter in 1:length(ListOfDatasets)){

    dataset <- ListOfDatasets[[counter]]

    if (verbose == TRUE) {
      cat("#################### GENERAL INFO ####################")
      cat("\n")
      cat("\n")
      print( paste( "Time Series",counter,"of", length(ListOfDatasets),". It starts on", index(head(dataset)[1]), "and ends on", index(tail(dataset)[6]) ), quote=FALSE, row.names = FALSE )
      cat("\n")
      print(median( diff( index(dataset) )))
      cat("\n")

      cat("Summary statistics:")
      print(summary(dataset))
      cat("\n")

      #Detect outliers
      cat("#################### SCANNING FOR OUTLIERS ####################")
      cat("\n")
      print(chisq.out.test(as.numeric(dataset)))
      cat("\n")
    }

    # detect the most probable time step
    timestep <- median( as.numeric(diff( index(dataset) )) )

    if (returnGapsInfo==TRUE & verbose == TRUE){
      cat("#################### SCANNING FOR GAPS IN THE RECORDS ####################")
      cat("\n")
      cat("\n")

      lengthTotal <-length(dataset)
      lengthMissing <- length(which(is.na(dataset)))
      percentageMissing <- round( (1 - ((lengthTotal-lengthMissing) / lengthTotal))*100, 2)

      cat("Percentage of missing values")
      print( paste(percentageMissing,"%", sep="") )
      if ( any(diff(index(dataset))>timestep) ) {
        cat("Missing values have been detected at the following starting times:")
        cat("\n")
        cat("\n")
        cat("2*timestep gap:")
        cat("\n")
        print( paste(index( dataset[which(diff(index(dataset))>timestep & diff(index(dataset))<=2*timestep)] ), ",", sep="" ), quote = FALSE)
        cat("\n")
        cat("3*timestep gap:")
        cat("\n")
        print( index( dataset[which(diff(index(dataset))>2*timestep & diff(index(dataset))<=3*timestep)] ), quote = FALSE)
        cat("\n")
        cat("Longer gaps:")
        cat("\n")
        print( index( dataset[which(diff(index(dataset))>3*timestep)] ), quote = FALSE)
        cat("\n")
      }else{
        cat( "No gaps have been detected!")
        cat("\n")
        cat("\n")
      }
    }

    if (returnTimeInfo==TRUE & verbose == TRUE){
      cat("#################### SCANNING FOR DIFFERENT RECORDING TIMES ####################")
      cat("\n")
      cat("\n")
      regularTS <- NA
      timestep <- as.numeric(timestep)
      # Check if this is a regular hourly time series
      if (attributes(median( diff( index(dataset) )))$units =="hours" & median( as.numeric(diff( index(dataset) )))==1) {
        timestep <- 60
        print("This is a regular hourly time series")
        regularTS <- "YES"
      }

      # check if there are missing values at the beginning of the dataset
      if (length( sort( unique(as.POSIXlt(index(dataset)[1:60/timestep])$min) ) ) < 60/timestep ) {
        cat("Missing values at the beginning of the dataset, please remove them and try again.")
      }else{

        tempdataset <- dataset
        RecordingTimes <- sort( unique(as.POSIXlt(index(tempdataset)[1:60/timestep])$min) )

        if (length(RecordingTimes)==4){

          RecordingChange <- which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[4])[1]
          if (!is.na(RecordingChange)) {
            cat( "This is an irregular time series, as recording times change throughout the dataset:" )
            cat("\n")
            cat("\n")
            cat( paste("From ",index(tempdataset)[1]," to ",index(tempdataset)[RecordingChange-1]," values are recorded every hour at the following minutes: ",
                       RecordingTimes[1], " ", RecordingTimes[2], " ",
                       RecordingTimes[3], " ", RecordingTimes[4], ".",sep="" ) )
            cat("\n")
            while ( any( as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[4] ) ){

              tempdataset <- tempdataset[RecordingChange:length(tempdataset)]
              RecordingTimes <- sort( unique(as.POSIXlt(index(tempdataset)[1:60/timestep])$min) )
              RecordingChange <- ifelse(is.na(which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[4])[1]),
                                        length(tempdataset),
                                        which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[4])[1])
              cat( paste("From ",index(tempdataset)[1]," to ",index(tempdataset)[RecordingChange-1]," values are recorded every hour at the following minutes: ",
                         RecordingTimes[1], " ", RecordingTimes[2], " ",
                         RecordingTimes[3], " ", RecordingTimes[4], ".",sep="" ) )
              cat("\n")
            }
          }else{
            cat( "This is a regular time series" )
            regularTS <- "YES"
            cat("\n")
          }
        }

        if (length(RecordingTimes)==6){
          RecordingChange <- which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[4] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[5] &
                                     as.POSIXlt(index(tempdataset))$min != RecordingTimes[6])[1]
          if (!is.na(RecordingChange)) {
            cat( "This is an irregular time series, as recording times change throughout the dataset:" )
            cat("\n")
            cat("\n")

            cat( paste("From ",index(tempdataset)[1]," to ",index(tempdataset)[RecordingChange-1]," values are recorded every hour at the following minutes: ",
                       RecordingTimes[1], " ", RecordingTimes[2], " ",
                       RecordingTimes[3], " ", RecordingTimes[4], " ",
                       RecordingTimes[5], " ", RecordingTimes[6], ".",sep="" ) )
            cat("\n")
            while ( any( as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[4] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[5] &
                           as.POSIXlt(index(tempdataset))$min != RecordingTimes[6] ) ){

              tempdataset <- tempdataset[RecordingChange:length(tempdataset)]
              RecordingTimes <- sort( unique(as.POSIXlt(index(tempdataset)[1:60/timestep])$min) )
              RecordingChange <- ifelse(is.na(which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[4] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[5] &
                                                      as.POSIXlt(index(tempdataset))$min != RecordingTimes[6])[1]),
                                        length(tempdataset),
                                        which(as.POSIXlt(index(tempdataset))$min != RecordingTimes[1] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[2] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[3] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[4] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[5] &
                                                as.POSIXlt(index(tempdataset))$min != RecordingTimes[6])[1])
              cat( paste("From ",index(tempdataset)[1]," to ",index(tempdataset)[RecordingChange-1]," values are recorded every hour at the following minutes: ",
                         RecordingTimes[1], " ", RecordingTimes[2], " ",
                         RecordingTimes[3], " ", RecordingTimes[4], " ",
                         RecordingTimes[5], " ", RecordingTimes[6], ".",sep="" ) )
              cat("\n")
            }
          }else{
            cat( "This is a regular time series" )
            cat("\n")
          }
        }

      }
      cat("\n")
    }

    if (returnNegInfo==TRUE & verbose == TRUE){

      cat("#################### SCANNING FOR NEGATIVE VALUES ####################")
      cat("\n")
      cat("\n")
      percentageNegative <- 0
      if (any(dataset<0,rm.na=TRUE)) {
        percentageNegative <- round( 100 * (1 - (lengthTotal-length(which(dataset<0)))/lengthTotal ), 2)
        if (percentageNegative > 0){
          startignore <- head(dataset[which(dataset<0)])[1]
          endignore <- tail(dataset[which(dataset<0)])[6]
          print(paste(percentageNegative,"% negative values detected between",index(startignore),"and", index(endignore)))
        }
      }else{
        cat( "No negative values have been detected!" )
        cat("\n")
        cat("\n")
      }
    }

  }

  summaryRecord <- list("Start"=index(head(dataset)[1]),
                     "End"=index(tail(dataset)[6]),
                     "Regular"=regularTS,
                     "percentageMissing"=percentageMissing,
                     "percentageNegative"=percentageNegative)

  return(summaryRecord)

}

#' Scans time series in a named list
#'
#' @param myTS list of time series, it can be zoo, xts, etc
#'
#' @return print useful information
#'
#' @examples
#' # myTS( dataset=Q[[8]])
#'

scanTimeSeries <- function(myTS){

  dimnamesTS <- attributes(myTS)[2]$dimnames[[2]]

  if ( is.null(dimnamesTS) ){

    results <- ScanTS(myTS)

  }else{

    lTS <- length(dimnamesTS)
    results <- list()

    for (i in 1:lTS){

      eval(parse(text=paste("myTS$",dimnamesTS[i],sep="")))

      results[1]

    }

  }

  return(results)

}

