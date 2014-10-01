#' GLUE uncertainty analysis
#'
#' @param DATA This is a data.frame containing the observed time series (zoo objects).
#'             It is structured into three columns:
#'             "P" containing precipitation,
#'             "E" containing evapo-transpiration and
#'             "Q" containing streamflow discharge.
#' @param deltim time step in days
#' @param models list of FUSE model structures to use
#' @param parameters table containing 24 parameters (columns) and a number of rows equal to the number of samples
#' @param threshold this is typically a number in the range [0.1,0.9]
#'
#' @return a plot showing uncertainties
#'
#' @author Claudia Vitolo
#'
#' @examples
#' # x <- glue(DATA, deltim, models, parameters, threshold)

glue <- function(DATA, deltim, models, parameters, threshold){

  warmup <- round(dim(DATA)[1]*0.1,0)
  lengthResults <- dim(DATA)[1] - warmup

  tableResults <- matrix(NA,ncol=3+lengthResults,nrow=0)

  Qobs <- DATA[(warmup+1):dim(DATA)[1],"Q"]

  for (mid in models) {

    print(paste(mid,"out of",length(models)))

    for ( pid in 1:dim(parameters)[1] ) {

      temp <- RunFUSE(DATA, parameters[pid,], deltim, mid)

      Qsim <- temp[(warmup+1):dim(DATA)[1]]

      NS <- EF(Qobs,Qsim)

      if (NS > threshold){
        newrow <- c(mid, pid, NS, t(Qsim))
        tableResults <- rbind(tableResults, newrow)
      }

    }
  }

  #Uncertainty
  NSeff <- tableResults[,3]
  Qsim  <- tableResults[,4:dim(tableResults)[2]]

  # Define Lower and Upper quantile for the prediction bounds.
  # Here we take the 0.05 and 0.95 quantiles resulting in 90% prediction limits.
  lower <- 0.05
  upper <- 0.95

  #Normalise the efficiencies so that they sum up to 1:
  eff <- NSeff - threshold
  eff <- eff/sum(eff)

  # Now we calculate the quantiles for each timestep
  limits <- apply( Qsim,2,"wtd.quantile", weights =  eff, probs = c(lower,upper), normwt=T)

  DateTime <- index(DATA[(warmup+1):dim(DATA)[1],])

  Qobs <- zoo(DATA[(warmup+1):dim(DATA)[1],"Q"],order.by=DateTime)
  Llimit <- zoo(as.numeric(limits[1,]),order.by=DateTime)
  Ulimit <- zoo(as.numeric(limits[2,]),order.by=DateTime)

  df2plot <- data.frame(DateTime,Qobs,Llimit,Ulimit)

  df <- melt(df2plot ,  id = 'DateTime', variable_name = 'variable')

  # plot on same grid, each series colored differently --
  # good if the series have same scale
  #   p <- ggplot(df, aes(DateTime,value)) +
  #     geom_line(aes(colour = variable)) +
  #     ggtitle(paste("Model ID =",mid))
  p <- ggplot(data=df2plot, aes(x = DateTime)) +
    geom_ribbon(aes(ymin=Llimit,ymax=Ulimit), alpha=0.35) +
    geom_line(aes(y=Qobs, color="blue")) +
    xlab(" ") +
    scale_y_continuous("Streamflow Discharge") +
    scale_colour_manual(values = "lightskyblue4") +
    scale_fill_manual("", values = c("coral4", "darkkhaki")) +
    ggtitle(paste("Model ID =",paste(models))) +
    theme(legend.title = element_blank())


  # or plot on different plots
  # ggplot(df, aes(DateTime,value)) + geom_line() + facet_grid(variable ~ .)

  return(list("t"=tableResults,"p"=p))

}
