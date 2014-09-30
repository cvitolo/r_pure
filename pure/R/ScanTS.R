#' Scan time series given in input and check for un-even recording times,
#' negative values and missing values.
#'
#' @param dataset Time series, it can be zoo, xts, etc
#' @param returnNegInfo bolean, default value is TRUE
#' @param returnTimeInfo bolean, default value is TRUE
#' @param returnMissingsInfo bolean, default value is TRUE
#'
#' @param verbose bolean, default value is TRUE
#'
#' @return print useful information and suggestions
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # ScanTS( dataset )
#'

ScanTS <- function(dataset,
                   returnNegInfo=TRUE,
                   returnTimeInfo=TRUE,
                   returnMissingsInfo=TRUE,
                   verbose = FALSE){

  # GENERAL INFO
  NEG  <- FALSE
  GAP  <- FALSE
  IRR  <- FALSE
  MIS  <- FALSE

  lengthTotal <- length(dataset)

  # detect the most probable time step
  timestep <- median( as.numeric(diff( index(dataset) )) )
  myUnits <- attributes(median( diff( index(dataset) )))$units

  if ( verbose == TRUE ){

    cat("#### GENERAL INFO #####################################################")
    cat("\n")
    cat("\n")

    print( paste( "This time series starts on",
                  index(head(dataset)[1]),
                  "and ends on",
                  index(tail(dataset)[6]) ),
           quote=FALSE, row.names = FALSE )
    cat("\n")

    print( paste("The most probable time step is",
                 timestep,
                 myUnits),
           quote=FALSE )
    cat("\n")

    cat("Summary statistics:")
    cat("\n")
    cat("\n")
    print(summary(dataset))
    cat("\n")

  }

  # NEGATIVE VALUES

  percentageNegative <- 0

  if ( any(dataset < 0) ) {

    lNeg <- length( which( dataset < 0 ) )
    percentageNegative <- 1 - ((lengthTotal - lNeg)/lengthTotal)

    if (percentageNegative > 0){

      startignore <- head(dataset[which(dataset<0)])[1]
      endignore <- tail(dataset[which(dataset<0)])[lNeg]

    }

    if ( returnNegInfo == TRUE ) NEG <- TRUE

  }

  if ( returnNegInfo==TRUE & verbose == TRUE ){

    cat("#### NEGATIVE VALUES ############################################")
    cat("\n")
    cat("\n")

    if (percentageNegative > 0){

      print(paste(percentageNegative,
                  "% negative values detected between ",
                  index(startignore),
                  " and ",
                  index(endignore),
                  sep=""))

    }else{

      cat( "No negative values have been detected!" )
      cat("\n")
      cat("\n")

    }

  }

  # GAPS IN RECORDS AND TIMESTEP IRREGULARITIES

  if ( any( diff(index(dataset)) > timestep ) ) {
    IRR <- TRUE

    tempdataset <- dataset
    tempIndex0 <- as.POSIXlt(index(tempdataset)[1:60/timestep])$min
    RecordingTimes <- sort( unique(tempIndex0) )
    tempIndex <- as.POSIXlt(index(tempdataset))$min

  }

  if ( returnTimeInfo==TRUE & verbose == TRUE ){
    cat("#### GAPS IN THE RECORDS ############################################")
    cat("\n")
    cat("\n")

    if ( any( diff(index(dataset)) > timestep ) ) {

      GAP <- TRUE

      cat("Gaps have been detected at the following starting times:")
      cat("\n")
      cat("\n")
      tstep <- diff(index(dataset))

      cat("2*timestep gap:")
      cat("\n")
      checkCondition <- which( tstep > timestep & tstep <= 2*timestep )
      print( paste(index( dataset[checkCondition] ), ",",sep=""),
             quote = FALSE )
      cat("\n")

      cat("3*timestep gap:")
      cat("\n")
      checkCondition <- which( tstep > 2*timestep & tstep <= 3*timestep )
      print( paste(index( dataset[checkCondition] ), ",",sep=""),
             quote = FALSE )
      cat("\n")

      cat("Longer gaps:")
      cat("\n")
      checkCondition <- which( tstep > 3*timestep )
      print( paste(index( dataset[checkCondition] ), ",",sep=""),
             quote = FALSE )
      cat("\n")

    }else{

      cat( "No gaps have been detected!")
      cat("\n")
      cat("\n")

    }

    if ( length(RecordingTimes) == 4 ){

      RecordingChange <- which(tempIndex != RecordingTimes[1] &
                                 tempIndex != RecordingTimes[2] &
                                 tempIndex != RecordingTimes[3] &
                                 tempIndex != RecordingTimes[4])[1]

      if (!is.na(RecordingChange)) {
        cat( "Recording times change throughout the dataset:" )
        cat("\n")
        cat("\n")
        cat( paste("From",index(tempdataset)[1],
                   "to",index(tempdataset)[RecordingChange-1]))
        cat("\n")
        cat( paste("values are recorded every hour at the following minutes:",
                   " ", RecordingTimes[1], " ", RecordingTimes[2], " ",
                   RecordingTimes[3], " ", RecordingTimes[4]) )
        cat("\n")
        cat("\n")
        while ( any( tempIndex != RecordingTimes[1] &
                       tempIndex != RecordingTimes[2] &
                       tempIndex != RecordingTimes[3] &
                       tempIndex != RecordingTimes[4] ) ){

          tempdataset <- tempdataset[RecordingChange:length(tempdataset)]

          tempIndex0 <- as.POSIXlt(index(tempdataset)[1:60/timestep])$min
          RecordingTimes <- sort( unique(tempIndex0) )
          tempIndex <- as.POSIXlt(index(tempdataset))$min

          x1 <- which(tempIndex != RecordingTimes[1] &
                        tempIndex != RecordingTimes[2] &
                        tempIndex != RecordingTimes[3] &
                        tempIndex != RecordingTimes[4])[1]

          RecordingChange <- ifelse(is.na(x1), length(tempdataset), x1)

          cat( paste("From",index(tempdataset)[1],
                     "to",index(tempdataset)[RecordingChange-1]))
          cat("\n")
          cat( paste("values are recorded every hour at the following minutes:",
                     " ", RecordingTimes[1], " ", RecordingTimes[2], " ",
                     RecordingTimes[3], " ", RecordingTimes[4], ".",sep="" ) )
          cat("\n")
          cat("\n")
        }

      }

    }

    if ( length(RecordingTimes) == 6 ){

      RecordingChange <- which(tempIndex != RecordingTimes[1] &
                                 tempIndex != RecordingTimes[2] &
                                 tempIndex != RecordingTimes[3] &
                                 tempIndex != RecordingTimes[4] &
                                 tempIndex != RecordingTimes[5] &
                                 tempIndex != RecordingTimes[6])[1]

      if (!is.na(RecordingChange)) {

        cat( "Recording times change throughout the dataset:" )
        cat("\n")
        cat("\n")

        cat( paste("From", index(tempdataset)[1],
                   "to", index(tempdataset)[RecordingChange-1]))
        cat("\n")
        cat( paste("values are recorded every hour at the following minutes:",
                   " ", RecordingTimes[1], " ", RecordingTimes[2], " ",
                   RecordingTimes[3], " ", RecordingTimes[4], " ",
                   RecordingTimes[5], " ", RecordingTimes[6], ".",sep="" ) )
        cat("\n")
        cat("\n")
        while ( any( tempIndex != RecordingTimes[1] &
                       tempIndex != RecordingTimes[2] &
                       tempIndex != RecordingTimes[3] &
                       tempIndex != RecordingTimes[4] &
                       tempIndex != RecordingTimes[5] &
                       tempIndex != RecordingTimes[6] ) ){

          tempdataset <- tempdataset[RecordingChange:length(tempdataset)]

          tempIndex0 <- as.POSIXlt(index(tempdataset)[1:60/timestep])$min
          RecordingTimes <- sort( unique(tempIndex0) )
          tempIndex <- as.POSIXlt(index(tempdataset))$min

          x1 <- which(tempIndex != RecordingTimes[1] &
                        tempIndex != RecordingTimes[2] &
                        tempIndex != RecordingTimes[3] &
                        tempIndex != RecordingTimes[4] &
                        tempIndex != RecordingTimes[5] &
                        tempIndex != RecordingTimes[6])[1]

          RecordingChange <- ifelse(is.na(x1), length(tempdataset), x1)

          cat( paste("From",index(tempdataset)[1],
                     "to",index(tempdataset)[RecordingChange-1]))
          cat("\n")
          cat( paste("values are recorded every hour at the following minutes:",
                     " ", RecordingTimes[1], " ", RecordingTimes[2], " ",
                     RecordingTimes[3], " ", RecordingTimes[4], " ",
                     RecordingTimes[5], " ", RecordingTimes[6], ".",sep="" ) )
          cat("\n")
          cat("\n")
        }

      }

    }

  }

  # MISSING VALUES

  lengthMissing <- length(which(is.na(dataset)))
  percentageMissing <- (1 - ((lengthTotal-lengthMissing) / lengthTotal))

  if ( percentageMissing > 0) {

    MIS <- TRUE

  }

  if ( returnMissingsInfo == TRUE & verbose == TRUE ){

    cat("#### MISSING VALUES #################################################")
    cat("\n")
    cat("\n")

    if (percentageNegative > 0){

      print( paste("Missing values", percentageMissing, "%"), quote=FALSE)

    }else{

      cat( "No missing values have been detected!" )
      cat("\n")
      cat("\n")

    }

  }

  message( "#### Recommendations: ###########################################" )


  message( paste("Do I need to correct negative values?",
                 ifelse(NEG,"YES","NO") ) )
  if ( NEG == TRUE ) {
    message("Use CorrectNeg() to remove negative values,")
    message("this will generate missing values!")
  }
  cat("\n")

  message( paste("Do I need to correct recording times to make them regular?",
                 ifelse(IRR,"YES","NO") ) )
  if ( GAP == TRUE & IRR == TRUE ) {
    message("Use FillGaps() then Irr2Reg() to get a regular time series")
  }
  if ( GAP == TRUE & IRR == FALSE ) {
    message("Use FillGaps() to add time steps for the detected gaps, ")
    message("this will generate missing values!")
  }
  if ( GAP == FALSE & IRR == TRUE ) {
    message("Use Irr2Reg() to get a regular time series")
  }
  cat("\n")

  message( paste("Do I need to infill missing values?",
                 ifelse(any(GAP,MIS),"YES","NO") ) )
  if ( any(GAP,MIS) == TRUE | all(GAP,MIS) == TRUE ) {
    message("Use FillGaps() to infill missing values")
  }
  cat("\n")

  message("Run the function ScanTS() with option verbose = TRUE for more info.")

}
