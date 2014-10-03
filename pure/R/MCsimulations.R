#' This function runs Monte Carlo simulations for a given number of realisations.
#'
#' @param DATA this is the data frame containing the observations
#' @param deltim this is the time step (in days)
#' @param warmup this is the warmup period
#' @param parameters table containing the parameter sets, see \code{GeneratePsetsFUSE}
#' @param ModelList this is the reduced list of model structures (total number: 240)
#' @param SimulationFolder path to the folder where the function saves the results
#' @param MPIs list of functions describing the Model performance Indices
#' @param verbose if TRUE, it prints details of the current simulation
#'
#' @return 2 arrays: R (discharges) and S (indices).
#'
#' @examples
#' # MCsimulations(DATA,deltim,warmup,parameters,ModelList,SimulationFolder,MPIs)
#'

MCsimulations <- function(DATA,deltim,warmup,parameters,ModelList,
                          SimulationFolder, MPIs,
                          verbose=TRUE){

  # library(fuse)
  DATA <- coredata(DATA)

  for (i in 1:dim(ModelList)[1]){ ###!!!

    mid <- ModelList[i,"mid"]
    indices <- data.frame(matrix(NA, nrow=dim(parameters)[1], ncol=length(MPIs) ))
    discharges <- data.frame(matrix(NA, nrow=dim(parameters)[1],ncol=dim(DATA)[1]-warmup))

    for (pid in 1:dim(parameters)[1]){

      if (verbose==TRUE) print(paste("Current model: ",mid," - Current run(p. set): ",pid))

      ParameterSet <- as.list(parameters[pid,])

      q_routed <- RunFUSE(DATA, ParameterSet, deltim, mid)

      x <- data.frame( Po = DATA[(warmup + 1):dim(DATA)[1],"P"],
                       Qo = DATA[(warmup + 1):dim(DATA)[1],"Q"],
                       Qs = q_routed[(warmup + 1):dim(DATA)[1]] )

      y <- lapply(MPIs, function(f) sapply(list(x), function(d) f(d) ) )

      indices[pid,] <- as.numeric(as.character(y))

      discharges[pid,] <- q_routed[(warmup + 1):dim(DATA)[1]]

    }

    names(indices) <- names(MPIs)

    save(indices, discharges, file = paste(SimulationFolder,"/MID_",mid,".Rdata",sep=""))

  }

}

