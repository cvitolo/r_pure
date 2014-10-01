#Emulation examples!

#The framework we're looking at involves
#1) A reality y
#2) An observation z on that reality made with error
#so z = y + observationError
#3) A computer simulator f(x) with inputs x
#4) A "best inputs" x^* such that f(x^*) is as close to
#y as it's possible to get from the simulator:
#y = f(x^*) + discrepancy
#5) Discrepancy can be split into internalDisc and externalDisc
#but that isn't too important right now

#Instead of doing examples using FUSE we'll start with some simpler examples

#Let us suppose our simulator is y(x) = exp(x) + 3*x

simulator <- function(x) exp(x) + 3*x

#Let us suppose that:
#The reality we have observed is z=15
data <- 15

#There is no error on the observation z
observationError <- 0

#Our simulator perfectly predicts reality
internalDisc <- 0
externalDisc <- 0

#This is the simplest case

#Suppose we think that bounds for x are [0,3]
lowerbounds <- 0
upperbounds <- 3

#Note that x^* here is about 2.15, let's see how close
#we can get to this.

#We want to try to emulate the output of this function
emulatorQuantityFunctions <- list(function(x) x)
#In more complicated situations, such as FUSE, we can't directly emulate
#the whole function so will be emulating a collection of summaries (e.g. maximum, average)
#Each element of emulatorQuantityFunctions will correspond to one of these summaries,
#so fx will always be a list of functions that will be applied to the model
#output to compute the summaries
#In this case, emulatorQuantityFunctions has only one element and that is the identity function

relevantParams <- 1
#When running FUSE, depending on mid some of the parameters have no effect.
#It makes no sense to vary them when creating the model runs
#relevantParams will in general be a vector that specifies which of the parameters
#are relevant. In this case there is only one parameter that is relevant, so we
#use a length 1 vector with value 1.

#Suppose we have run the function at x=0,1,2,3
#The functions tend to like matrices with each row corresponding to all the parameters
#for a given run. Since we only have one parameter per run, we have one column in the matrix
params <- matrix(0:3, ncol=1)
#You can just give the functions a vector but then they tend
#to assume that you have one run with n parameters rather than
#n runs with 1 parameter
modelRuns <- as.vector(simulator(params))
#we want modelRuns in vector form though!


#The first step for emulation is to work out a good fit.
#An emulator is a simple approximation for the model in question
#The emulator will be of the form f(x) = \sum \beta_i g_i(x) + u(x)
#u is the local variation
#\sum \beta_i g_i(x) gives the mean trend
#we specify the g_i as a list of functions
#Let's try a linear fit as the simplest case.
meanFunction <- list(meanFunctionLinear)
#In general we can build up as complicated a mean function as we like

#Let's see how good the fit is
stepwiseEmulatorSelection(params, modelRuns, emulatorQuantityFunctions[[1]], meanFunction)
#You'd do this once for each emulatorQuantityFunction element;
#you'd want a different mean function for each emulatorQuantityFunction element
#In this case the interesting number is [2,1] which is the adjusted R^2 of this fit
#which isn't bad but let's see if we can do better.
#stepwiseEmulatorSelection has other powers that we'll see later

meanFunction <- list(meanFunctionQuadratic, meanFunctionLinear)
stepwiseEmulatorSelection(params, modelRuns, emulatorQuantityFunctions[[1]], meanFunction)

#OK, that looks like a good fit! We'll use that
#typically this step can be rather harder but this examples was chosen rather nicely.
#Now we're in a position to emulate.
#We'll use a gaussian covariance function for u
covarianceFunction <- list(covarianceGaussian)
#Again, when there are multiple elements to emulatorQuantityFunction
#we'll have multiple covariance functions
#We'll use a correlation parameter of 1
delta <- 1
#but we might change that later.

#Let's have a look at the emulator for
newPoints <- matrix(seq(0, 3, 0.1), ncol=1)

#The next function takes a new collection of parameters newPoints
#and gives the emulator's expectation E(simulator(newPoints))
#and the variance Var(simulator(newPoints))
emulatorTest <- emulatorExpAndVarEmpirical(newPoints, params, relevantParams, modelRuns, emulatorQuantityFunctions,
                                           list(meanFunction), list(delta), covFunction=covarianceFunction)
#here we have lists for meanFunction delta because in general there needs to be
#one of them for every element of covarianceFunction
#this will give a warning because relevantParams "ought" to be a matrix but it doesn't matter here
#The output of emulatorExpAndVar is a 3-dimensional array.
dim(emulatorTest)
#The first dimension corresponds to the different emulatorQuantityFunctions (only one here).
#The second dimension corresponds to each of the new points
#The third dimension corresponds to expectation (1) and variance(2)
#So emulatorTest[1,65,2] gives Var(emulatorQuantityFunctions[[1]](simulator(newPoints[65,])))
#That is, the variance of the emulator for emulator quantity 1 at the 65th parameter choice
#in newPoints.

#now let's look at our emulator
plot(newPoints, emulatorTest[1,,1], type='l')
#this is the prediction (expected value) for the simulator's output at each point
lines(newPoints, emulatorTest[1,,1] + 3*sqrt(abs(emulatorTest[1,,2])), col='blue')
lines(newPoints, emulatorTest[1,,1] - 3*sqrt(abs(emulatorTest[1,,2])), col='blue')
#3 standard deviations (abs is because occasionally when asking the emulator
#for the variance at points we've already sampled it gives thing like -10^-14 instead
#of the correct 0)
lines(newPoints, simulator(newPoints), col='red')
#actual simulator output

#So it looks pretty good, although for x near 3 we're hitting the edge of 3sds; this could
#be a sign that our correlation distance is a bit dodgy

#Now the question is: what possible values might give us a value of 15?
#We can use the plots we've created to basically work out what the possible range is,
#but in more complicated situations we won't have such a simple visualization and pattern
secondWave <- seq(0, 3, 0.01)
#a large collection of points
secondWave <- secondWave[-c(1,101,201,301)]
#remove the points we've already sampled; we know because of perfect information and perfect model
#that they aren't good enough.
secondWave <- matrix(secondWave, ncol=1)
#as matrix again

historyMatch(secondWave, params, modelRuns, rep(data, 297), 1, 
             emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta))
#History match returns all the points in secondWave whose emulator's expected value
#was within 3 standard deviations of observed z
#The standard deviation is calculated by summing the variances of internalDisc, externalDisc
#and observationError (all zero here) along with the emulator variance at each point.

#Let's look a bit more closely at the internal workings of historyMatch.
#We can get more information that we might be able to visualize nicely.
#The first call is emulatorDiagnostics.
diagnostics <- emulatorDiagnostics(params, modelRuns,
                    secondWave, 1,
                    emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                    internalDisc=0, 
                    externalDisc=0, 
                    observationError=0, sds=3)
diagnostics$lower
diagnostics$upper
#These are the emulator expectations +/- 3*variance for each point in secondWave.
#We can use these in two ways:
#1) if we have the actual simulator runs for each element, we can test whether or not
#the simulator runs lie within these bounds. If most of them do, then everything is OK
#If many do not, then our emulator is over-confident
#If we have very many points and all of them are within the range then our emulator may
#be being too conservative

checkEmulatorDiagnostics(diagnostics, rep(data, nrow(secondWave)), emulatorQuantityFunctions)
#This has three parts.
#The first is a matrix where i,j is true if observed value of emulatorQuantityFunctions[[i]]
#is within the bounds for secondWave[j,].
#The second part is a vector that tells you how many of emulatorQuantityFunctions were within
#the bounds for each parameter choice (here there is only one element of emulatorQuantityFunctions
#so the values can only be zero or 1)
#The third part tells you which rows of secondWave were within the bounds for all emulatorQuantityFunctions
#historyMatch then reports just the third part.

#There is a second use for checkEmulatorDiagnostics
checkEmulatorDiagnostics(diagnostics, apply(secondWave, 1, simulator), emulatorQuantityFunctions)
#This is now checking whethe rthe simulator's value at each point is within the bounds.
#This is useful for emulator validation (hence the name "diagnostics")

#Important distinction: when doing HISTORY MATCHING you need to include
#observationError, internalDisc, externalDisc because we are matching the emulator to REALITY
#When doing DIAGNOSTICS you must not include these three elements (apart from nugget effect---see later)
#because diagnostics is about matching the emulator to the SIMULATOR.

#Normally we would now take the points that survived history matching
#and construct a new emulator using them.
#But in this case we can see that the correct answer lies within say [2.1, 2.2]
#We can do a bunch of runs on this subset
paramsSecondWave <- seq(2.1, 2.2, 0.01)
modelRunsSecondWave <- simulator(paramsSecondWave)
paramsSecondWave <- matrix(paramsSecondWave, ncol=1)

plot(paramsSecondWave, modelRunsSecondWave, type='l')
#looks like a linear fit will do now!
meanFunctionSecondWave <- list(meanFunctionLinear)

stepwiseEmulatorSelection(paramsSecondWave, modelRunsSecondWave, emulatorQuantityFunctions[[1]], meanFunctionSecondWave)

#same covariance function
covarianceFunctionSecondWave <- covarianceFunction
#but different delta (we're looking in a smaller region, so want to allow more local variation)
deltaSecondWave <- 1/20

#Create a large number of new points.
newPoints <- matrix(seq(2.1, 2.2, 0.001), ncol=1)

#First we history match using the first wave
thirdWaveRetainedFirstWave <- historyMatch(newPoints, params, modelRuns, rep(data, nrow(newPoints)), 1, 
                                       emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta))

#and now second wave
thirdWaveRetainedSecondWave <- historyMatch(thirdWaveRetainedFirstWave, paramsSecondWave, modelRunsSecondWave, rep(data, nrow(thirdWaveRetainedFirstWave)), 1, 
                                           emulatorQuantityFunctions, list(meanFunctionSecondWave), covarianceFunctionSecondWave, list(deltaSecondWave))
#We get an alert about negative variances, but don't worry, it's just an irrelevance here

thirdWaveRetainedSecondWave
#Empty! All our points were judged implausible.
#Let's have a look and see what's happening

emulatorTest <- emulatorExpAndVarEmpirical(newPoints, paramsSecondWave, relevantParams, modelRunsSecondWave,
                                           emulatorQuantityFunctions,
                                           list(meanFunctionSecondWave), list(deltaSecondWave), 
                                           covFunction=covarianceFunctionSecondWave)
#now let's look at our emulator
plot(newPoints, emulatorTest[1,,1], type='l')
#this is the prediction (expected value) for the simulator's output at each point
lines(newPoints, emulatorTest[1,,1] + 3*sqrt(abs(emulatorTest[1,,2])), col='blue')
lines(newPoints, emulatorTest[1,,1] - 3*sqrt(abs(emulatorTest[1,,2])), col='blue')
lines(newPoints, simulator(newPoints), col='red')
#Essentially our emulator is predicting near-perfectly, so unless we input the actual root
#of simulator(x) = 15, we will get 0 non-implausible points.
#Let's throw in the correct value to see what we get
trueRoot <- matrix(c(2.146987,2.146988,2.146989),ncol=1)
#The middle value is the true root
historyMatch(trueRoot, paramsSecondWave, modelRunsSecondWave, rep(data,3), 1, 
                                            emulatorQuantityFunctions, list(meanFunctionSecondWave), covarianceFunctionSecondWave, list(deltaSecondWave))
#Yes! We successfully judged the real answer non-implausible. Success!

#right, so what we have here is a very elaborate method to find the roots of equations...
#But it really becomes useful when we have some uncertainties.
#Now let's see what happens if we add some discrepancies. Suppose that we believe there is a
#normal error with 95% of readings within 5% of reality y. A variance of
observationError <- (3/8)^2
#models this.

#Let's also suppose that we have estimated an internal discrepancy of
internalDisc <- 0.2

#We'll keep the external discrepancy as zero (there isn't actually any distinction between
#the errors in the process anyway at this point)

#There's no difference in fitting the emulator: the simulator is
#still the same, it's just the link between the output and reality that has changed
#So all that changes is the history match step

historyMatch(secondWave, params, modelRuns, rep(data, 297), 1, 
             emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
             internalDisc=internalDisc, observationError=observationError)
#Observe that there are more non-implausible points now

paramsSecondWave <- seq(1.98, 2.28, 0.01)
modelRunsSecondWave <- simulator(paramsSecondWave)
paramsSecondWave <- matrix(paramsSecondWave, ncol=1)

plot(paramsSecondWave, modelRunsSecondWave, type='l')
#looks like a linear fit will do now!
meanFunctionSecondWave <- list(meanFunctionLinear)

stepwiseEmulatorSelection(paramsSecondWave, modelRunsSecondWave, emulatorQuantityFunctions[[1]], meanFunctionSecondWave)

#same covariance function
covarianceFunctionSecondWave <- covarianceFunction
#but different delta
deltaSecondWave <- 1/30

paramsThirdWave <- matrix(seq(1.98, 2.28, 0.001),ncol=1)
paramsThirdWave <- paramsThirdWave[-21,,drop=FALSE]

thirdWaveRetainedFirstWave <- historyMatch(paramsThirdWave, params, modelRuns, rep(data, nrow(paramsThirdWave)), 1, 
                                           emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                                           internalDisc=internalDisc, observationError=observationError)

thirdWaveRetainedFirstWave

thirdWaveRetainedSecondWave <- historyMatch(thirdWaveRetainedFirstWave, paramsSecondWave, modelRunsSecondWave, rep(data, nrow(thirdWaveRetainedFirstWave)), 1, 
             emulatorQuantityFunctions, list(meanFunctionSecondWave), covarianceFunctionSecondWave, list(deltaSecondWave),
             internalDisc=internalDisc, observationError=observationError)

thirdWaveRetainedSecondWave

#Observe that now the second wave adds no useful information
#The first wave tells us everything we can know given the discrepancies.

#Note that we can do similar sorts of plots for history matching as we did for emulator performance
#This is for the first wave:
newPoints <- matrix(seq(0, 3, 0.1), ncol=1)

emulatorTest <- emulatorExpAndVarEmpirical(newPoints, params, relevantParams, modelRuns, emulatorQuantityFunctions,
                                           list(meanFunction), list(delta), covFunction=covarianceFunction)
plot(newPoints, emulatorTest[1,,1], type='l')
lines(newPoints, emulatorTest[1,,1] + 3*sqrt(abs(emulatorTest[1,,2])+internalDisc+observationError), col='blue')
lines(newPoints, emulatorTest[1,,1] - 3*sqrt(abs(emulatorTest[1,,2])+internalDisc+observationError), col='blue')
lines(newPoints, rep(15, nrow(newPoints)), lty=2, col='red')

#####################
#END OF FIRST EXAMPLE
#####################

#Now let's think about a more complicated simulator with multiple inputs but only one output
simulator <- function(x) {exp(x[1]) + 3*x[1] - x[1]*x[2] + log(x[2]) - x[3]*10^(-6)}

#for x in [0,3], y in [2,4], z in [0,5].

#Let us suppose that:
#The reality we have observed is z=15
data <- 15

#There is no error on the observation z
observationError <- 0

#Our simulator perfectly predicts reality
internalDisc <- 0
externalDisc <- 0

#Bounds on the parameters
lowerbounds <- c(0,2,0)
upperbounds <- c(3,4,5)
#These become important this time.

#We want to try to emulate the output of this function
emulatorQuantityFunctions <- list(function(x) x)

relevantParams <- c(1,1,1)
#again, all three parameters are required to run the function.

#This time we have to work a bit harder to get a set of parameters to run the model at.
#Let's run it at 50 points.
params <- modelrunHypercube(50, relevantParams, lowerbounds, upperbounds)
params

modelRuns <- apply(params, 1, simulator)
emulatorQuantityFunctions <- list(function(x) x)

#When we have parameters with different ranges, it becomes useful to
#transform all of them onto a standard range, for instance 0 to 1
#This is because the emulator depends on the distances between different parameter choices
#and if we don't standardise the ranges then parameters with very big ranges
#will be considered distant even if they are close and parameters with very small ranges
#will be considered close even if they are distant. I have seen this destroy an emulator
#in horrible ways.
paramsStandardised <- standardiseParams(params, relevantParams, lowerbounds, upperbounds, c(0,1))
#the main pitfall from now on is remembering when you need your parameters standardised and
#when you don't!

meanFunction <- list(meanFunctionLinear)
#Let's see how good the fit is
stepwiseEmulatorSelection(paramsStandardised, modelRuns, emulatorQuantityFunctions[[1]], meanFunction)
#The output is larger than before. Again, [2,1] has the R^2 for the
#initial fit, but now there is another column next to it:
#[1,2] suggests the index of the variable that is doing the least in the fit
#and [2,2] is the R^2 with this variable removed
#The third column is an artefact of bad coding.

#So this suggests that we might want to remove variable 3 from consideration
#Getting rid of "inactive" parameters like x_3 can make things work better.
#Note that we can't completely remove it, because later in the process
#it might become active again. So, when creating new runs for later waves
#we do need to vary x_3.
activeParams <- c(1,1,0)

#Whenever we remove parameters we have to consider how much extra variance
#this might be introducing into our emulator.
#I haven't got a function to do this for any given function yet
#so we need to define one here.

estimateNuggetEffect <- function(params, activeParams, repetitions)
{
  #First of all insert activeParams into the full set of parameters
  #needed to run FUSE
  fullParams <- params
  #if fullParams isn't a matrix, make it one
  if(is.matrix(fullParams)==FALSE)
  {
    fullParams <- matrix(fullParams, nrow=1)
  }
  
  #each row of this matrix gives a combination of indices for
  #params, mid, fracstate0.
  
  #Create the inactiveParams
  inactiveParams <- relevantParams - activeParams
  inactiveParamsHypercube <- modelrunHypercube(repetitions, inactiveParams, lowerbounds, upperbounds) #hypercube
  SDs <- rep(0, nrow(fullParams))
  #Loop through the rows of inputIndices
  for(i in 1:nrow(fullParams))
  {
    thisParam <- fullParams[i,]
    
    thisParam <- insertParameterValues(inactiveParams, 
                                       inactiveParamsHypercube,
                                       baseParams=thisParam)
    
    output <- apply(thisParam, 1, simulator)
    SDs[i] <- sd(output)
  }
  return(SDs)
}

#What this function does: when given a collection of parameters params,
#for each row in params (set of parameter choices) it varies the inactive parameters
#repetitions times, and calculates the standard deviation
#So you get a standard deviation for each row of params.
#If this standard deviation is small and is similar for each row, then
#1) the choice of removing the inactive parameters seems to make sense
#2) the standard deviation shoudl be added to the internal discrepancy
nug <- estimateNuggetEffect(paramsStandardised, activeParams, 100)
internalDisc <- internalDisc + mean(nug)^2

#Unsurprisingly the example was constructed so that removing this parameter
#does essentially nothing.

#Now, our original linear fit wasn't particularly good.
#How about quadratic with interactions?
meanFunction <- list(meanFunctionQuadratic, meanFunctionLinearInteraction, meanFunctionLinear)

#This looks better, and parameter x_3 is still inactive.
#So, on we go!

covarianceFunction <- list(covarianceGaussian)
delta <- 1/3 #This is usually a good starting value for the correlation distance on [0,1]
#but of course we might tweak this later.

#We've got 50 runs, so we can do a test using 40 of them to see how well we're doing
diagnostic <- emulatorDiagnostics(paramsStandardised[1:40,], modelRuns[1:40],
                                  paramsStandardised[41:50,],
                                matrix(c(1,1,0), nrow=1),
                                emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                                internalDisc=internalDisc)

checkEmulatorDiagnostics(diagnostic, modelRuns[41:50], emulatorQuantityFunctions)

#I didn't get a great performance here, only 5 out of 10 within the bounds
#Chances are our correlation distance is a bit big

#If you're not particularly interested in the practical details of fitting
#a good emulator, feel free to skip this bit until END SKIP

delta <- 1/4
diagnostic <- emulatorDiagnostics(paramsStandardised[1:40,], modelRuns[1:40],
                                  paramsStandardised[41:50,],
                                  matrix(c(1,1,0), nrow=1),
                                  emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                                  internalDisc=internalDisc)

checkEmulatorDiagnostics(diagnostic, modelRuns[41:50], emulatorQuantityFunctions)

#You see that fitting a good emulator isn't an automatic procedure and you have
#to play around a lot to get a good answer...

#Of course, when we have a fast model we can check our emulator rather more crudely...
hypcube <- cbind(modelrunHypercube(1000, activeParams, lowerbounds, upperbounds), rep(1, 1000))
hypcubeStandardised <- standardiseParams(hypcube, relevantParams, lowerbounds, upperbounds, c(0,1))
diagnostic <- emulatorDiagnostics(paramsStandardised, modelRuns,
                                  hypcubeStandardised,
                                  matrix(c(1,1,0), nrow=1),
                                  emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                                  internalDisc=internalDisc)

testRuns <- apply(hypcube, 1, simulator)

sum(checkEmulatorDiagnostics(diagnostic, testRuns, emulatorQuantityFunctions)$numWithinBounds)
#I got 850 with delta=1/3
#930 with delta=1/4
#let's go with 1/4
delta <- 1/4
diagnostic <- emulatorDiagnostics(paramsStandardised, modelRuns,
                                  hypcubeStandardised,
                                  matrix(c(1,1,0), nrow=1),
                                  emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),
                                  internalDisc=internalDisc)

testRuns <- apply(hypcube, 1, simulator)

sum(checkEmulatorDiagnostics(diagnostic, testRuns, emulatorQuantityFunctions)$numWithinBounds)
#Do the bad runs have anything in common?
badRuns <- (1:1000)[checkEmulatorDiagnostics(diagnostic, testRuns, emulatorQuantityFunctions)$numWithinBounds ==0]
hypcube[badRuns,]
#Yes: the emulator performs poorly at the extreme of x[1].
#This may be something to worry about:
testRuns[badRuns]
#some of these are close to reality.
#So some plausible points might be removed.
#we could consider countering this by adding a few extra runs in the areas we're doing badly in.
#But let's go with making delta 1/5

#END SKIP

delta <- 1/5

#Now let's do a history match
secondWaveParams <- modelrunHypercube(1000, c(1,1,1), lowerbounds, upperbounds)
secondWaveParamsStandardised <- standardiseParams(secondWaveParams, relevantParams, lowerbounds, upperbounds, c(0,1))
#a large collection of points
#Note that we include x[3] in here: even though it is inactive in the first emulator
#it might be active in the second wave emulator and so we need a spread of these points

secondWaveRetained <- historyMatch(secondWaveParamsStandardised, paramsStandardised, modelRuns, rep(data, 1000), matrix(activeParams, nrow=1), 
             emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta), internalDisc=internalDisc)
secondWaveRetained
#I got 8 surviving

modelRunsSecondWave <- apply(unstandardiseParams(secondWaveRetained, c(1,1,1), lowerbounds, upperbounds), 1, simulator)
modelRunsSecondWave
#looks OK!

plot(secondWaveParams[,1:2])
points(unstandardiseParams(secondWaveRetained, c(1,1,1), lowerbounds, upperbounds), col='red')
#A potential worry is that we're removing too much here.
#In cases such as this where we're a bit worried about our emulator we can
#accept points up to say 4 standard deviations from reality rather than the
#usual 3.
#We won't try that here though.

#Now, one way to proceed would be to use these remaining points to fit a new emulator.
#But there are only 8.
#We could also make our hypercube bigger and bigger until we get a good number
#Alternatively, note 
unstandardiseParams(secondWaveRetained, c(1,1,1), lowerbounds, upperbounds)
#suggests say 2.4 to 2.8 for x[1]
#We want to be fairly conservative though, so let's say
lowerboundsSecondWave <- c(2.2,2,0)
upperboundsSecondWave <- c(3, 4, 5)

#Let's make a big hypercube in this region
secondWaveParams <- modelrunHypercube(1000, c(1,1,1), lowerboundsSecondWave, upperboundsSecondWave)

#We want to run the first emulator on all of these, so we need to standardise them in the INITIAL way
secondWaveStandardised <- standardiseParams(secondWaveParams, c(1,1,1), lowerbounds, upperbounds)

secondWaveRetained <- historyMatch(secondWaveStandardised, paramsStandardised, modelRuns, rep(data, 1000), matrix(activeParams, nrow=1), 
                                   emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta),internalDisc=internalDisc)

plot(secondWaveParams[,1:2])
points(unstandardiseParams(secondWaveRetained, c(1,1,1), lowerbounds, upperbounds)[,1:2], col='red')

#This looks like we've roughly worked out the trend
secondWaveRuns <- apply(unstandardiseParams(secondWaveRetained, c(1,1,1), lowerbounds, upperbounds), 1, simulator)
secondWaveRuns
#Looks good

stepwiseEmulatorSelection(secondWaveRetained, secondWaveRuns, emulatorQuantityFunctions[[1]], meanFunction)
#Looks like 3 is still having no real effect.

#We could now go on to do everything again but I won't bother
#Essentially the next wave destroys anything not within 0.002 of 15.

#Now let us imagine
observationError <- (3/4)^2
internalDisc <- internalDisc + 0.2

#Let's compare history matching with no discrepancies
secondWaveParams <- modelrunHypercube(1000, c(1,1,1), lowerbounds, upperbounds)
secondWaveParamsStandardised <- standardiseParams(secondWaveParams, relevantParams, lowerbounds, upperbounds, c(0,1))
secondWaveRetained <- historyMatch(secondWaveParamsStandardised, paramsStandardised, modelRuns, rep(data, 1000), matrix(activeParams, nrow=1), 
                                   emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta), internalDisc=internalDisc-0.2)
plot(secondWaveParamsStandardised[,1:2])
points(secondWaveRetained[,1:2], col='red')

#with discrepancies.
secondWaveRetainedDisc <- historyMatch(secondWaveParamsStandardised, paramsStandardised, modelRuns, rep(data, 1000), matrix(activeParams, nrow=1), 
                                   emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta), internalDisc=internalDisc,
                                       observationError=observationError)
plot(secondWaveParamsStandardised[,1:2])
points(secondWaveRetainedDisc[,1:2], col='blue')

#Observe how many more points we retained when adding this discrepancy!

#We might be able to get away with a linear fit on this reduced region
meanFunctionSecondWave <- list(meanFunctionLinear)
secondWaveRunsDisc <- apply(unstandardiseParams(secondWaveRetainedDisc, c(1,1,1), lowerbounds, upperbounds),
                            1, simulator)

stepwiseEmulatorSelection(secondWaveRetainedDisc, secondWaveRunsDisc, emulatorQuantityFunctions[[1]], meanFunctionSecondWave)
#Looks good!


#Try a third wave
thirdWaveParams <- modelrunHypercube(5000, c(1,1,1), lowerbounds, upperbounds)
thirdWaveParamsStandardised <- standardiseParams(thirdWaveParams, relevantParams, lowerbounds, upperbounds, c(0,1))

thirdWaveRetainedFirstWave <- historyMatch(thirdWaveParamsStandardised, paramsStandardised, modelRuns, rep(data, 5000), matrix(activeParams, nrow=1), 
                                       emulatorQuantityFunctions, list(meanFunction), covarianceFunction, list(delta), internalDisc=internalDisc,
                                       observationError=observationError)

thirdWaveRetainedSecondWave <- historyMatch(thirdWaveRetainedFirstWave, secondWaveRetainedDisc, secondWaveRunsDisc, rep(data, nrow(thirdWaveRetainedFirstWave)), matrix(activeParams, nrow=1), 
                                           emulatorQuantityFunctions, list(meanFunctionSecondWave), covarianceFunction, list(delta), internalDisc=internalDisc,
                                           observationError=observationError)

plot(thirdWaveParamsStandardised[,1:2])
points(thirdWaveRetainedFirstWave[,1:2], col='blue')
points(thirdWaveRetainedSecondWave[,1:2], col='red')
#That region looks like it's all down to the discrepancy now and we won't be able to do any better.
#Indeed, if you plot the function in this reduced region you find that the function is essentially a hyperplane and
#so our linear fit will be essentially the best thing we can ever manage.