#' This function runs Monte Carlo simulations for a given number of realisations.
#'
#' @param DATA this is the data frame containing the observations
#' @param deltim this is the time step (in days)
#' @param warmup this is the warmup period
#' @param parameters table containing the parameter sets, see \code{GeneratePsetsFUSE}
#' @param ModelList this is the reduced list of model structures (total number: 240)
#' @param SimulationFolder path to the folder where the function saves the results
#' @param verbose if TRUE, it prints details of the current simulation
#'
#' @return 2 arrays: R (discharges) and S (indices).
#'
#' @examples
#' # MCsimulations(DATA,deltim,warmup,parameters,ModelList,SimulationFolder)
#'

MCsimulations <- function(DATA,deltim,warmup,parameters,ModelList,
                          SimulationFolder,
                          verbose=TRUE){

  # library(fuse)
  DATA <- coredata(DATA)

  for (i in 1:dim(ModelList)[1]){ ###!!!

    mid <- ModelList[i,"mid"]
    indices <- data.frame(matrix(NA, nrow=dim(parameters)[1],ncol=5))
    discharges <- data.frame(matrix(NA, nrow=dim(parameters)[1],ncol=dim(DATA)[1]-warmup))

    for (pid in 1:dim(parameters)[1]){

      if (verbose==TRUE) print(paste("Current model: ",mid," - Current run(p. set): ",pid))

      ParameterSet <- as.list(parameters[pid,])

      q_routed <- RunFUSE(DATA, ParameterSet, deltim, mid)

      indices[pid,] <- PerformanceIndices(Po=DATA[(warmup + 1):dim(DATA)[1],"P"],
                                          Qo=DATA[(warmup + 1):dim(DATA)[1],"Q"],
                                          Qs=q_routed[(warmup + 1):dim(DATA)[1]])

      discharges[pid,] <- q_routed[(warmup + 1):dim(DATA)[1]]

    }

    names(indices) <- c("LAGTIME","MAE","NSHF","NSLF","RR")

    save(indices, discharges, file = paste(SimulationFolder,"/MID_",mid,".Rdata",sep=""))

  }

}

