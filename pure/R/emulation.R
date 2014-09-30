#A collection of functions to facilitate history matching


#The first bunch of functions are going to be related to transformations
#between the actual parameter values and a standardised range for each one
#boundaries for parameter values and the extent of the standardised range should
#be user configurable.

#So to use this function: params should be a vector or matrix of parameters
#relevantParams should be a boolean vector as follows:
#there will be n parameters that are necessary to specify (22 for FUSE).
#In a given moment we might only be varying m<=n parameters and keeping the
#rest fixed.
#We specify a standard ordering for our parameters (for FUSE I use the table in
#Clarke et al that describes the parameters as the ordering)
#relevantParams is then a vector with m 1s, and the position of the 1s determines
#which m of the n parameters we're actually looking at.
#So if the first 1 was in position 3, that would mean that the first element/column
#of params corresponds to parameter 3 in the model
#upperbounds and lowerbounds are the upper and lower bounds for the n parameters
#in the case of FUSE have used the bounds from Clarke et al, but made the time delay
#upper bound lower because in the Pontbren catchment the higher values for time delay
#are always poor.
#standardisedBounds should be a vector of 2 elements giving the lower and upper bounds
#for the standardised parameters.
standardiseParams <- function(params,relevantParams,
                              lowerbounds=c(50,100,0.05,0.05,0.05,0.05,0.01,1,1,1,0.05,0.01,0.001,1,0.001,0.001,0.001,0.05,0.001,5,2,0.01),
                              upperbounds=c(5000,10000,0.95,0.95,0.95,0.95,1000,20,250,5,0.95,1000,10000,10,0.25,0.25,0.25,0.95,3,10,5,2),
                              standardisedBounds=c(0,1))
{
  #starting checks
  paramsUsed <- params
  if(is.matrix(paramsUsed)==FALSE)
  {
    paramsUsed <- matrix(paramsUsed, nrow=1)
    warning("params was not a matrix; made it one with 1 row, check this makes sense")
  }

  if(sum(relevantParams)!=ncol(paramsUsed))
  {
    stop("number of relevantParams does not coincide with columns of params")
  }

  if(length(relevantParams)!=length(lowerbounds))
  {
    stop("number of parameters in relevantParams not the same as in lowerbounds")
  }

  if(length(relevantParams)!=length(upperbounds))
  {
    stop("number of parameters in relevantParams not the same as in upperbounds")
  }

  if(length(standardisedBounds)!=2)
  {
    stop("standardisedBounds does not have length 2")
  }

  upper <- upperbounds[relevantParams==1]
  lower <- lowerbounds[relevantParams==1]
  b <- -(-standardisedBounds[1]+(standardisedBounds[2]-standardisedBounds[1])*lower/(upper-lower))
  a <- (standardisedBounds[2]-standardisedBounds[1])/(upper-lower)
  ans <- matrix(ncol=ncol(paramsUsed), nrow=nrow(paramsUsed))
  for(i in 1:nrow(paramsUsed))
  {
    ans[i,] <- a*paramsUsed[i,] + b
  }
  return(ans)
}

unstandardiseParams <- function(params,relevantParams,
                                lowerbounds=c(50,100,0.05,0.05,0.05,0.05,0.01,1,1,1,0.05,0.01,0.001,1,0.001,0.001,0.001,0.05,0.001,5,2,0.01),
                                upperbounds=c(5000,10000,0.95,0.95,0.95,0.95,1000,20,250,5,0.95,1000,10000,10,0.25,0.25,0.25,0.95,3,10,5,2),
                                standardisedBounds=c(0,1))
{
  paramsUsed <- params
  if(is.matrix(paramsUsed)==FALSE)
  {
    paramsUsed <- matrix(paramsUsed, nrow=1)
    warning("params was not a matrix; made it one with 1 row, check this makes sense")
  }

  if(sum(relevantParams)!=ncol(paramsUsed))
  {
    stop("number of relevantParams does not coincide with columns of params")
  }

  if(length(relevantParams)!=length(lowerbounds))
  {
    stop("number of parameters in relevantParams not the same as in lowerbounds")
  }

  if(length(relevantParams)!=length(upperbounds))
  {
    stop("number of parameters in relevantParams not the same as in upperbounds")
  }

  if(length(standardisedBounds)!=2)
  {
    stop("standardisedBounds does not have length 2")
  }

  upper <- upperbounds[relevantParams==1]
  lower <- lowerbounds[relevantParams==1]
  b <- -(-standardisedBounds[1]+(standardisedBounds[2]-standardisedBounds[1])*lower/(upper-lower))
  a <- (standardisedBounds[2]-standardisedBounds[1])/(upper-lower)
  ans <- matrix(ncol=ncol(paramsUsed), nrow=nrow(paramsUsed))
  for(i in 1:nrow(paramsUsed))
  {
    ans[i,]<-(paramsUsed[i,]-b)/a
  }
  return(ans)
}

#Now we use these to create a standardised hypercube to run the model or do history mathcing with

modelrunHypercube <- function(repetitions, relevantParams,
                                lowerbounds = c(50,100,0.05,0.05,0.05,0.05,0.01,1,1,1,0.05,0.01,0.001,1,0.001,0.001,0.001,0.05,0.001,5,2,0.01),
                                upperbounds = c(5000,10000,0.95,0.95,0.95,0.95,1000,20,250,5,0.95,1000,10000,10,0.25,0.25,0.25,0.95,3,10,5,2))
{
  #params: order as in Clark et al.
  #upper and lower bounds
  #NOTE: timedelay range much smaller than in Clark because most of Clark range seems to give silly answers
  #lowerparams <- c(50,100,0.05,0.05,0.05,0.05,0.01,1,1,1,0.05,0.01,0.001,1,0.001,0.001,0.001,0.05,0.001,5,2,0.01)
  #upperparams <- c(5000,10000,0.95,0.95,0.95,0.95,1000,20,250,5,0.95,1000,10000,10,0.25,0.25,0.25,0.95,3,10,5,5)

  lowerParams <- lowerbounds[relevantParams==1]
  upperParams <- upperbounds[relevantParams==1]

  numparams <- length(lowerParams)
  hypercube<-maximinLHS(repetitions, numparams)
  #unstandardise
  cube <- unstandardiseParams(hypercube, relevantParams, lowerbounds, upperbounds, standardisedBounds=c(0,1))
  return(cube)
}

#The next function looks at how good a particular mean function is for fitting an emulator,
#and also suggests which parameters could be treated as irrelevant.
#params is a matrix or vector of model parameters at which we've run the model
#modelRuns is a collection (of the appropriate type--see later) of model runs, one for each
#row of param.
#emulatorQuantityFunction is a function that will be applied to modelRuns to get whatever
#emulator quantity is needed. It must take as input modelRuns, so the sort of object
#modelRuns should be is determined by this function.
#For instance for FUSE, modelRuns could be a zoo object with a time series for each row of params,
#an emulatorQuantityFunction could be a function that takes a zoo object and returns the maximum
#of each time series in that object.
#If you have already applied such a function to modelRuns, then for this function set modelRuns to
#be the output of this function and emulatorQuantityFunction to function(x) {x}.
#Finally, meanFunction is a list of functions that together create the mean function desired for the emulator.
#For instance, it could be a list contained a meanFunctionLinear, meanFunctionQuadratic, meanFunctionInteractions.
#The function then fits param ~ emulatorQuantityFunction(modelRuns) according to this mean function, and
#determines the adjusted R^2.
#The next step is to consider removing each of the parameters in turn and run the fit again.
#The parameter that reduces the R^2 by the least is then removed, and this new R^2 recorded.
#This process continues until there are no more parameters left.
#The output is a matrix whose first row corresponds to the parameter index removed at that stage
#and the second row corresponds to the adjusted R^2 from this model.
#If the first R^2 isn't high, then the mean function is a bad fit
#We should then consider treating as irrelevant all parameters that could be removed without
#damaging the R^2 too much.
#It's important to note that we shouldn't jut take the output of this function in an
#automatic way and instead think about what it tells us before using its advice or ignoring it.
stepwiseEmulatorSelection <- function(params, modelRuns, emulatorQuantityFunction, meanFunction)
{
  #in standard computer model notation, we have simulator f that acts on inputs x (params)
  #and returns the quantity of interest.
  #In this case the quantity of interest is emulatorQuantityFunction(modelRuns), so let's
  #calculate that

  paramsUsed <- params
  #if params is not a matrix, make it one
  if(is.matrix(params)==FALSE)
  {
    paramsUsed <- matrix(params, nrow=1)
    warning("params was not a matrix, made it one but check this is OK")
  }

  fx <- do.call(emulatorQuantityFunction, args=list(modelRuns))

  #We start by creating the fit that we're interested in and recording the R^2.

  #columnRemoved is going to record which column we're considering removing next
  columnRemoved <- rep(0, ncol(paramsUsed))
  #Rsquared is going to be the Rsquared from removing that column
  Rsquared <- columnRemoved

  #with nothing removed, create the fit.
  #The fit is of the form \sum \beta_i meanFunction[[i]](x)
  #So we need to know the meanFunction at each x

  meanFunctionSimulated <- matrix(nrow=nrow(paramsUsed), ncol=0)
  for(i in 1:length(meanFunction))
  {
    #looping is the most straightforward way to do this
    meanFunctionSimulated <- cbind(meanFunctionSimulated, do.call(meanFunction[[i]], list(paramsUsed)))
  }

  #So now we've got the meanFunctionSimulated for this fx
  #Now we need to construct the formula object to feed into lm.
  mainFit <- lm(fx ~ ., as.data.frame(meanFunctionSimulated))

  #record the R squared!
  Rsquared[1] <- summary(mainFit)$adj.r.squared

  #We're going to slowly take away columns of params and
  #elements of fx so we'll need to have objects to record there

  paramsReduced <- paramsUsed
  fxReduced <- fx

  #What we're now going to do is to loop through each
  #parameter and see what removing it would do to the fit.
  #Then we choose the one with the lowest change in R^2 and remove it.

  counter <- 2

  while(counter < ncol(paramsUsed))
  {
    bestRSquared <- 0
    currentBest <- 0
    for(i in 1:ncol(paramsReduced))
    {
      #remove column
      paramsTemp <- paramsReduced[,-i, drop=FALSE]
      #we need to recalculate emulSim (we could be cleverer but I can't be bothered)
      meanFunctionSimulatedReduced <- matrix(nrow=nrow(paramsTemp), ncol=0)
      for(j in 1:length(meanFunction))
      {
        #looping is the most straightforward way to do this
        meanFunctionSimulatedReduced <- cbind(meanFunctionSimulatedReduced, do.call(meanFunction[[j]],
                                                                                    list(paramsTemp)))
      }
      reducedFit <- lm(fx ~ ., as.data.frame(meanFunctionSimulatedReduced))

      #Is this fit better than the best we've recorded?
      if(summary(reducedFit)$adj.r.squared > bestRSquared)
      {
        #record the R squared!
        bestRSquared <- summary(reducedFit)$adj.r.squared
        currentBest <- i
      }
    }
    Rsquared[counter] <- bestRSquared
    columnRemoved[counter]<-currentBest
    paramsReduced <- paramsReduced[,-currentBest, drop=FALSE]
    counter <- counter+1
  }

  #Now let's try to get columnRemoved into a more sensible format.

  #the first row of the output gives the variables removed at each stage
  #the second row gives the R squared for the quadratic fit on the
  #reduced model at this stage.
  return(matrix(c(columnRemoved, Rsquared), nrow=2, byrow=TRUE))

}

#Identifies the model structure associated with switch choices.
#In order: upper layer, lower layer, surface, percolation,
#evaporation, interflow, timedelay
identifyModelStructure <- function (s2,s3,s4,s5,s6,s7,s8, modlist)
{
  #find matrix row that gives required model structure
  #modlist is a data frame where each column represents a particular switch.
  #the value of a particular element is xy where x is the column's number
  #and y is the number of the switch. So for each s_i we create a number
  #i s_i (that is, i*10 + s_i). All the rows of the matrix with value other
  #that this can be removed. Doing this for all columns at the same time
  #should leave us with one row.
  reals2 <- as.numeric(paste(as.character(2), as.character(s2), sep=""))
  reals3 <- as.numeric(paste(as.character(3), as.character(s3), sep=""))
  reals4 <- as.numeric(paste(as.character(4), as.character(s4), sep=""))
  reals5 <- as.numeric(paste(as.character(5), as.character(s5), sep=""))
  reals6 <- as.numeric(paste(as.character(6), as.character(s6), sep=""))
  reals7 <- as.numeric(paste(as.character(7), as.character(s7), sep=""))
  reals8 <- as.numeric(paste(as.character(8), as.character(s8), sep=""))

  subset(modlist, arch1 == reals2 & arch2==reals3 & qsurf==reals4 & qperc==reals5 &esoil==reals6 & qintf==reals7 & q_tdh==reals8)[1,1]

}

insertParameterValues <- function(relevantParams, insertedValues,
                                  baseParams=c(81,2500,0.9,0.6,0.4,0.8,10,13,20,2,0.9,400,700,6,0.2,0.2,0.2,0.9,1.4,7.5,4,2.5))
{
  #Suppose we've got a base collection of parameters, and a set of values
  #for the relevantParams that we want to insert into the base collection.
  #In practice, relevantParams will be the set of parameters that actually
  #influences our chosen model structure. The simulator doesn't run if
  #it isn't given all the parameters, but the non-relevantParams don't have
  #any influence. Therefore when running the model for different parameter choices,
  #we can keep the non relevantParams the same. So, we start with a base collection
  #baseParams and replace all the relevantParams with values.
  outputRows <- 1
  if(is.matrix(insertedValues))
  {
    #We want a row for each row of insertedValues
    outputRows <- nrow(insertedValues)
  }

  answer <- matrix(rep(baseParams, outputRows), nrow=outputRows, byrow=TRUE)

  indices <- (1:length(relevantParams))[relevantParams==1]
  answer[,indices]<-insertedValues
  return(answer)
}

#OK, so now we've got a function for doing the stepwise selection for
#each emulation quantity.
#This is going to give us a logical matrix activeParamaters
#It'd be nice for convenience to have a simple function that takes
#a matrix of points and strips out the irrelevant parameter columns
#This function takes a MATRIX of points and a VECTOR of activeVariables
#and returns a MATRIX
#So if we want all the matrices we have to apply

removeInactiveParameters <- function(params, activeParameters)
{
  #Do some checks

  paramsMatrix <- params
  if(is.matrix(paramsMatrix)==FALSE)
  {
    #let's hope it's a vector
    paramsMatrix <- matrix(points, nrow=1)
    warning("params wasn't a matrix, so I made it one. If it wasn't a vector to start with, things might go wrong")
  }


  if(ncol(paramsMatrix) != length(activeParameters))
  {
    stop("length of logical activeVariables vector differs from number of parameters in points")
  }

  return(paramsMatrix[,(1:length(activeParameters))[activeParameters==TRUE], drop=FALSE])
}

#Suppose we've got a collection of relevantParams and of those
#a subcollection of activeParams.
#We'd like to estimate how much the relevantParams that are not
#activeParams influence the emulatorQuantityFunctions
#We estimate this by: 1) construct a collection of combinations of
#activeParams, mids, and fracstates0s; 2) for each element of this
#collection, vary the non-activeParams; 3) for each variation, run the model
#4) for each element of the collection, calculate the sd of
#emulatorQuantityFunction for all the variations; 5) report the results.
estimateNuggetEffect <- function(params, relevantParams, activeParams,
                                 mid, fracstate0, P, E, deltim,
                                 emulatorQuantityFunction, repetitions,
                                 baseParams=c(81,2500,0.9,0.6,0.4,0.8,10,13,20,2,0.9,400,700,6,0.2,0.2,0.2,0.9,1.4,7.5,4,2.5))
{
  #P and E should be time series or vectors
  #fractstate0, mid can be vectors or scalars
  #params should be a matrix or a vector
  #relevantParams should be a logical vector
  #activeParams should be a logical vector
  #ncol(params) should be sum(activeParams)
  #sum(relevantParams) <= sum(activeParams)
  #repetitions is how many different choices of non-activeParams
  #we use.

  #First of all insert activeParams into the full set of parameters
  #needed to run FUSE
  fullParams <- insertParameterValues(activeParams, params, baseParams=baseParams)
  #if fullParams isn't a matrix, make it one
  if(is.matrix(fullParams)==FALSE)
  {
    fullParams <- matrix(fullParams, nrow=1)
  }

  #Now create all the possibilities that we're going to run things at
  inputIndices <- expand.grid(1:nrow(fullParams), 1:length(mid),
                              1:length(fracstate0))

  #each row of this matrix gives a combination of indices for
  #params, mid, fracstate0.

  #Create the inactiveParams
  inactiveParams <- relevantParams - activeParams
  inactiveParamsHypercube <- TODO #hypercube
  SDs <- rep(0, nrow(inputIndices))
  #Loop through the rows of inputIndices
  for(i in 1:nrow(inputIndices))
  {
    thisParam <- fullParams[inputIndices[i,1],]
    thisMid <- mid[inputIndices[i,2]]
    thisFracstate0 <- fracstate0[inputIndices[i,3]]

    thisParam <- insertParameterValues(inactiveParams,
                                       inactiveParamsHypercube,
                                       baseParams=thisParam)

    output <- RUNFUSE #TODO
    SDs[i] <- APPLYFUNCTION #TODO
  }
  return(SDs)
}

#So far we've been using activeParams as a vector the same length as relevantParams.
#When we're doing emulation and history matching, we don't need the non-relevantParams
#at all.
#The following convenience function turns an activeParams of length relevantParams
#to an activeParams of length sum(relevantParams)
reallyActiveParams <- function(relevantParams, activeParams)
{
  return(activeParams[relevantParams==1])
}

#Right, so now we have the tools to create a hypercube, run it through
#FUSE, try some fits, identify inactive parameters, estimate the
#nugget effects of those parameters.
#Next step is to create the emulatorExpAndVar functions,

emulatorExpAndVarEmpirical <- function(newPoints, oldPoints, activeParameters, modelRuns, fx, meanFunction, delta, outputOrder=c(3,2,1), fxRunsCalculated=FALSE, covFunction)
{
  #We are emulating k quantities that are computing by the functions
  #in the k-list fx.
  #We have run the model at oldPoints, and observed output modelRuns
  #Emulator quantity fx[[i]] will be emulated by meanFunction[[i]]
  #with correlation distance delta[i].
  #active parameters for quantity fx[[i]] are given by activeParameters[i,]
  #We calculate expectation and variance of fx for each point in
  #newPoints

  #if fxRunsCalculated is TRUE then fx will have already been applied
  #to modelRuns in which case modelRuns is an m-vector
  #if it's false, then modelRuns must be the output of the model
  #that fx can apply to
  #Set it to TRUE when modelRuns is something really big and hard to work
  #with and you don't want to keep it around,
  #or if fx is expensive

  #The output will be an array where one dimension is emulation quantity
  #one dimension corresponds to expectation or variance
  #one dimension corresponds to newPoints
  #the order of dimension is given by outputOrder
  #c(1,2,3) would give exp/var, newPoints, fx
  #default gives fx, newPoints, exp/var

  #Covariance function is given by covFunction; this needs to have the following arguments:
  #x (matrix of points), y (matrix of points), betaVar (vector), delta (vector), sigmasq (scalar),
  #emul.x (matrix), emul.y (matrix), distanceMatrix (matrix)
  #We also need to change delta to be a list of vectors (one element for each emulation quantity)
  #Right now it doesn't do anything because we haven't finished checking everything
  #Actually let's make this a list of functions, because then we can have a different covariance
  #structure for each quantity if we want.

  #Initial checks

  #making things Matrices that might not be matrices
  newPointsMatrix <- newPoints
  oldPointsMatrix <- oldPoints
  activeParametersMatrix <- activeParameters
  if(is.matrix(newPointsMatrix)==FALSE)
  {
    newPointsMatrix <- matrix(newPointsMatrix, nrow=1)
    warning("newPoints wasn't a matrix; I have made it one")
  }

  if(is.matrix(oldPointsMatrix)==FALSE)
  {
    oldPointsMatrix <- matrix(newPointsMatrix, nrow=1)
    warning("oldPoints wasn't a matrix; I have made it one")
  }

  if(is.matrix(activeParametersMatrix)==FALSE)
  {
    activeParametersMatrix <- matrix(activeParametersMatrix, nrow=1)
    warning("activeParameters wasn't a matrix; I have made it one")
  }

  #Checking lengths of things
  if(ncol(oldPointsMatrix)!=ncol(newPointsMatrix))
  {
    print(ncol(oldPointsMatrix))
    print(ncol(newPointsMatrix))
    stop("oldPoints and newPoints have different number of columns")
  }

  if(ncol(oldPointsMatrix)!=ncol(activeParametersMatrix))
  {
    print(ncol(oldPointsMatrix))
    print(ncol(activeParametersMatrix))
    stop("oldPoints and activeParameters have different number of columns")
  }

  if(length(fx)!=nrow(activeParametersMatrix))
  {
    print(length(fx))
    print(nrow(activeParametersMatrix))
    stop("activeParameters has different number of rows than the length of fx")
  }

  deltaList <- delta
  if(is.list(deltaList)==FALSE)
  {
    warning("delta is not a list, making into one")
    deltaList <- list(deltaList)
  }
  #if(length(delta)==1 && length(fx)>1)
  #{
  #  deltaVector <- rep(delta, length(fx))
  #  warning("delta has length 1 and fx is longer, replicating delta to length(fx)")
  #}

  if(length(fx) != length(deltaList))
  {
    print(length(fx))
    print(length(deltaVector))
    stop("fx and delta have different lengths")
  }

  if(length(fx) != length(meanFunction))
  {
    print(length(fx))
    print(length(meanFunction))
    stop("fx and meanFunction have different lengths")
  }

  #Function proper starts here

  #Let's calculate fx (modelRuns) (that is, fxSim), assuming fxRunsCalculated=FALSE

  if(fxRunsCalculated==FALSE)
  {
    fxSim <- sapply(fx, do.call, args=list(modelRuns))
    #this will be a matrix with columns for each element of fx;
    #we want this the other way round
    fxSim <- t(fxSim)
  }else
  {
    #fx of Runs have already been calculated, so no need to do it
    fxSim <- modelRuns
  }

  #if fxSim isn't a matrix, better make it one
  if(is.matrix(fxSim)==FALSE)
  {
    warning("fxSim isn't a matrix; if you have more than one emulation quantity then something is wrong")
    fxSim <- matrix(fxSim, nrow=1)
  }


  #For each fxSim[i], run the subfunction
  #We want to take with us: all of oldPoints, all of newPoints,
  #delta[i], meanFunction[[i]],
  #activeParameters[i,].
  #Clearly this is a suitable setting for mapply.
  #To use mapply we need to split activeParameters into a list
  activeParametersList <- split(activeParametersMatrix,
                                1:nrow(activeParametersMatrix))
  fxSimList <- split(fxSim, 1:nrow(fxSim))

  expAndVar <- mapply(emulatorExpAndVarEmpiricalSingleFx, fxSim=fxSimList,
                      delta=deltaList,meanFunction=meanFunction,
                      activeParameters=activeParametersList, covFunction=covFunction,
                      MoreArgs=list(oldPoints=oldPointsMatrix,
                                    newPoints=newPointsMatrix),
                      SIMPLIFY="array")

  #this will be a 3 dimensional array.
  #The first dimension will be length 2
  #The second dimension will be n
  #The third dimension will be k
  #so [1,j,k] will be the expected value for point j for emulator quantity k
  #[2,j,k] will be the variance for point j for emulator quantity k

  #I'd prefer it so that the dimensions are switched to k, n, 2
  #but there's no reason to no allow the user to specify this
  expAndVar <- aperm(expAndVar, outputOrder)
  return(expAndVar)
}

emulatorExpAndVarEmpiricalSingleFx <- function(newPoints, oldPoints,
                                               activeParameters, fxSim, meanFunction,
                                               delta, covFunction)
{
  #Now we have matrices newPoints and oldPoints,
  #logical vector activeParameters,
  #an emulator quantity fxSim
  #list of functions meanFunction
  #delta is some quantity that feeds into covFunction
  #and we are in a position to do the main chunk of the work.

  numModelRuns <- nrow(oldPoints)

  #Remove irrelevant parameters
  newPointsReduced <- removeInactiveParameters(newPoints, activeParameters)
  oldPointsReduced <- removeInactiveParameters(oldPoints, activeParameters)

  #Let's calculate how many betas we have for this variable
  #We'll do that by applying all the meanFunctions to one element of x
  numBeta <- length(constructEmulatorCoefficients(oldPointsReduced[1,,drop=FALSE],
                                                  meanFunction))+1
  #+1 because the intercept is not contained in constructEmulatorCoefficients

  #set up the things we need to build the emulator

  #create the multipliers of the beta for all oldPoints
  emulSim <- t(apply(oldPointsReduced, 1, constructEmulatorCoefficients,
                     meanFunction=meanFunction))
  if(nrow(emulSim)==1)
  {
    #might need to transpose; apply might have gone bad
    emulSim <- t(emulSim)
  }

  #get the dimnames right!
  colnames(emulSim) <- colnames(constructEmulatorCoefficients(oldPointsReduced[1,], meanFunction))
  #So now we have a matrix of the coefficients, and we can do the
  #regression on it

  fit <- lm(fxSim ~ ., as.data.frame(emulSim))

  betaPrior <- fit$coefficients
  betaVar <- (summary(fit)$coefficients[,2])^2
  sigmasq <- (summary(fit)$sigma)^2

  #Now we have to find the
  #prior expected values of oldPoints

  #add a column of ones at the start!
  emulSimComplete <- cbind(rep(1, nrow(emulSim)), emulSim)

  #We need the distance matrix for oldPointsReduced
  oldPointsDistance <- as.matrix(dist(oldPointsReduced))
  dimnames(oldPointsDistance) <- NULL
  distanceMatrixSim <- oldPointsDistance

  #prior variance matrix
  #Ainv <- solve(covarianceFunction(oldPointsReduced, oldPointsReduced,
  #                                 betaVar, delta, sigmasq))
  Ainv <- ginv(do.call(covFunction, args=list(x=oldPointsReduced, y=oldPointsReduced,
                                               betaVar=betaVar, delta=delta, sigmasq=sigmasq,
                                               emul.x=emulSimComplete, emul.y=emulSimComplete,
                                               distanceMatrix=distanceMatrixSim)))


  expValueSim <- betaPrior %*% t(emulSimComplete)
  expValueSim <- expValueSim[1,]
  #expValueSim <-apply(t(betaPrior * t(emulSimComplete)), 1, sum)

  #We need the distance matrix in both exp and var, so we
  #may as well do it here.

  combinedPoints <- rbind(newPointsReduced,oldPointsReduced)
  combinedDistance <- as.matrix(dist(combinedPoints))
  dimnames(combinedDistance) <- NULL
  distanceMatrix <- combinedDistance[1:(nrow(newPointsReduced)),
                                     (1+nrow(newPointsReduced)):(nrow(newPointsReduced)+nrow(oldPointsReduced)),
                                     drop=FALSE]

  #We seem to have acquired pretty much everything we need now.
  #So let's apply adjustedExpectation and adjustedVariance
  #to each row of newPoints
  #We need to carry the ith row of the distance matrix with
  #the ith row of newPoints
  #so this is a job for mapply

  #Of course, we need distanceMatrix to be a list rather than a matrix.
  distanceList <- split(distanceMatrix, 1:nrow(distanceMatrix))

  #split newPoints
  newPointsList <- split(newPointsReduced, 1:nrow(newPointsReduced))

  xExp <- mapply(adjustedExpectation, newPoint=newPointsList, distanceVector=distanceList,
                 MoreArgs=list(oldPoints=oldPointsReduced, fxSim=fxSim, emulSim=emulSimComplete,
                               expValueSim=expValueSim, betaPrior=betaPrior,
                               betaVar=betaVar, Ainv=Ainv, delta=delta,
                               sigmasq=sigmasq, meanFunction=meanFunction, covFunction=covFunction))

  xVar <- mapply(adjustedVariance, newPoint=newPointsList, distanceVector=distanceList,
                 MoreArgs=list(oldPoints=oldPoints, emulSim=emulSimComplete, betaVar=betaVar,
                               Ainv=Ainv, delta=delta, sigmasq=sigmasq, meanFunction=meanFunction,
                               covFunction=covFunction))

  return(rbind(xExp, xVar))
}

constructEmulatorCoefficients <- function(point, meanFunction)
{
  #Here we have point being a single parameter VECTOR, and
  #meanFunction is a LIST of functions that'll build up the
  #emulator mean function.

  #the function determines the g_i(point)s, that will multiply
  #the beta coefficients.

  emulSim <- matrix(nrow=1, ncol=0)
  for(i in 1:length(meanFunction))
  {
    emulSim <- cbind(emulSim, do.call(meanFunction[[i]], list(point)))
  }
  return(emulSim)
}

adjustedExpectation <- function(newPoint, oldPoints, fxSim, emulSim,
                                expValueSim, betaPrior, betaVar,Ainv,
                                delta, sigmasq, meanFunction, distanceVector, covFunction)
{
  #newPoint is a vector, oldPoints is a matrix
  #meanFunction is a list of functions to create the coeffs of the betas in the mean function

  #Let's do that first
  emulNew <- constructEmulatorCoefficients(newPoint, meanFunction)
  #add the 1 and make a matrix
  emulNew <- matrix(c(1, emulNew), nrow=1)

  #   #We also need the distance matrix
  #   combinedPoints <- rbind(newPoint, oldPoints, deparse.level=0)
  #   combinedDistance <- as.matrix(dist(combinedPoints))
  #   dimnames(combinedDistance)<-NULL
  #   distanceMatrix <- combinedDistance[1, 2:(1+nrow(oldPoints)), drop=FALSE]

  #Turn the distance vector into a 1-row matrix
  #Not necessary, can use a vector
  #distanceMatrix <- matrix(distanceVector, nrow=1)

  covar <- do.call(covFunction, args=list(x=matrix(newPoint, nrow=1), y=oldPoints,
                                          betaVar=betaVar, delta=delta, sigmasq=sigmasq,
                                          emul.x=emulNew, emul.y=emulSim,
                                          distanceMatrix=distanceVector))

  adjustedExpectation <- sum(betaPrior *emulNew) +
    covar %*% Ainv %*% (fxSim-expValueSim)

  return(adjustedExpectation)
}


adjustedVariance <- function(newPoint, oldPoints, emulSim, betaVar,Ainv,
                             delta, sigmasq, meanFunction, distanceVector, covFunction)
{
  #newPoint is a vector, oldPoints is a matrix
  #meanFunction is a list of functions to create the coeffs of the betas in the mean function

  #Let's do that first
  #Note: we need this for both so we could do it outside and then add it to the mapply
  emulNew <- constructEmulatorCoefficients(newPoint, meanFunction)
  emulNew <- matrix(c(1, emulNew),nrow=1)
  #Turn the distance vector into a 1-row matrix
  #no need, vector works fine
  #distanceMatrix <- matrix(distanceVector, nrow=1)

  #assuming betaVar is a diagonal matirx, that is, no correlation
  covt <- do.call(covFunction, args=list(x=matrix(newPoint, nrow=1), y=oldPoints,
                                         betaVar=betaVar, delta=delta, sigmasq=sigmasq,
                                         emul.x=emulNew, emul.y=emulSim,
                                         distanceMatrix=distanceVector))

  adj <- sum(betaVar * emulNew^2) + sigmasq - covt%*%Ainv%*%t(covt)
  return(adj)
}




#Let's think about emulator diagnostics.
#Simplest idea is to have a set of runs to build the emulator
#and a set of runs to test against.

#So, this function will take emulatorParams, diagnosticParams,
#activeParams, internalDiscrepancy, externalDiscrepancy, observationError,
#sds, and the emulation building stuff:
#covarianceFunction, meanFunction, delta, modelRuns, fx
#discrepancies and errors should be zero for emulator diagnostics;
#it's just the history matching that calls this function too for which it matters.
emulatorDiagnostics <- function(emulatorParams, emulatorRuns,
                                diagnosticParams,
                                activeParams,
                                fx, meanFunction, covarianceFunction, delta,
                                internalDisc=rep(0, length(fx)),
                                externalDisc=rep(0, length(fx)),
                                observationError=rep(0, length(fx)), sds=3)
{
  #emulator expectation and variance for diagnostics
  emulateDiagnostics <- emulatorExpAndVarEmpirical(diagnosticParams, emulatorParams, activeParams, emulatorRuns, fx,
                                                   meanFunction, delta, outputOrder=c(3,2,1), fxRunsCalculated=FALSE,
                                                   covFunction=covarianceFunction)

  #This will be of the following format:
  #An array where element i,j,k
  #is k=1: expectation for fx[[i]] at parameter j
  #k=2: variance for fx[[i]] at parameter j

  #Next we work out the emulator expectation plus minus sds worth of
  #standard deviations
  diagnosticsExpectations <- emulateDiagnostics[,,1]
  diagnosticsVariances <- emulateDiagnostics[,,2]
  if(is.array(diagnosticsExpectations)==FALSE)
  {
    diagnosticsExpectations <- matrix(diagnosticsExpectations, nrow=1)
  }
  if(is.array(diagnosticsVariances)==FALSE)
  {
    diagnosticsVariances <- matrix(diagnosticsVariances, nrow=1)
  }


  variances <- diagnosticsVariances + matrix(rep(internalDisc + externalDisc + observationError,
                                                          ncol(diagnosticsExpectations)), nrow=1)

  if(min(variances) < 0)
  {
    warning("Alert: negative variances, could just be irrelevant computation error, please check")
    print(min(variances))
  }

  diagnosticLowerInterval <- diagnosticsExpectations - sds*sqrt(abs(variances))
  diagnosticUpperInterval <- diagnosticsExpectations + sds*sqrt(abs(variances))
  return(list(lower=diagnosticLowerInterval, upper=diagnosticUpperInterval))
}

#Then we can use diagnostic model runs to see how many and which ones
#are outside the boundaries.
checkEmulatorDiagnostics <- function(diagnostic, modelRuns, fx)
{
  fxDiagnostics <- sapply(fx, do.call, args=list(modelRuns))

  #This should form a matrix where the columns are the fx
  #and the rows are the fx(modelRuns)

  withinBounds <- diagnostic$lower <= t(fxDiagnostics) & diagnostic$upper >= t(fxDiagnostics)

  #This is a matrix where i,j is true if the emulator for quantity i
  #and run j is within the bounds.

  numWithinBounds <- apply(withinBounds, 2, sum)
  #This is a vector giving the number of quantities that were
  #within the bounds for each diagnostic run

  isWithinBounds <- (1:ncol(withinBounds))[numWithinBounds==nrow(withinBounds)]
  #This is the indices of the runs that were within the bounds for
  #all emulator quantities

  return(list(withinBounds=withinBounds, numWithinBounds=numWithinBounds,
              isWithinBounds=isWithinBounds))

}

#Now how about a history match procedure.
#Let's suppose we have a standardised hypercube
#and a collection of runs, etc.
#Then this procedure isn't too difficult.

historyMatch <- function(hypercube, emulatorParams, modelRuns, data,
                         activeParams, fx, meanFunction, covarianceFunction, delta,
                         internalDisc=rep(0, length(fx)),
                         externalDisc=rep(0, length(fx)),
                         observationError=rep(0, length(fx)), sds=3)
{
  # A simple application of the diagnostic functions will do the job
  emulatorPredictions <- emulatorDiagnostics(emulatorParams, modelRuns,
                                             hypercube,
                                             activeParams, fx, meanFunction,
                                             covarianceFunction, delta,
                                             internalDisc,
                                             externalDisc,
                                             observationError, sds)

  historyMatchResults <- checkEmulatorDiagnostics(emulatorPredictions, data, fx)
  return(hypercube[historyMatchResults$isWithinBounds,,drop=FALSE])
}

#This will return the portion of the hypercube that is not implausible
#for any f_i(x).
