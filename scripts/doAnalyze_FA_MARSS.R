
require(MASS)
require(MARSS)
require(ini)
# require(reshape2)
source("../projectSrc/utils/kalmanFilter/buildGoNogoAndLaserInputs.R")
source("../projectSrc/utils/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/fit_FA_MARSS.R")
source("../commonSrc/stats/kalmanFilter/create_FA_MARSS.R")

processAll <- function() {
    region <- "v1"
    shaft <- 1
    obsDim <- 7
    obsInputMemorySecs <- 0.0
    stateDim <- rep(3, times=length(obsInputMemorySecs))
    stateInputMemorySecs <- obsInputMemorySecs
    analysisStart <- 180 #sec
    analysisDur <- 60 # sec
    maxIter <- 20
    kfFunc <- "MARSSkfss"
    stateOffsetType <- "unconstrained"
    stateCovType <- "diagonal and unequal"
    obsOffsetType <- "unconstrained"
    obsCovType <- "diagonal and unequal"
    initialStateMeanType <- "unconstrained"
    initialStateCovType <- "diagonal and unequal"
    dataFilename <- "results/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_timeSeries.RData"
    resultsFilenamePattern <- "results/analysis_%s_shaft%d_obsDim%d_stateDim%02d_stateInputMemory%.02f_obsInputMemory%.02fstart%.02fsec_dur%.02fsec.RData"

    stopifnot(length(stateDim)==length(stateInputMemorySecs) & length(stateInputMemorySecs)==length(obsInputMemorySecs))
    #
    data <- get(load(dataFilename))
    sRate <- data$sRate
    binSizeSecs <- (data$breaks[2]-data$breaks[1])/sRate
    analysisSamples <- (analysisStart*sRate)+(1:(analysisDur*sRate))
    spikeCounts <- data[[sprintf("%sShaft%dSpikeCounts", region, shaft)]][1:obsDim,analysisSamples]
    goStim <- data$goStim[analysisSamples]
    nogoStim <- data$nogoStim[analysisSamples]
    laserStim <- data$laserStim[analysisSamples]
    spikeRates <- spikeCounts/binSizeSecs
    tSpikeRates <- sqrt(spikeRates)

    for(i in 1:length(stateDim)) {
        show(sprintf("Processing number of latents=%d, state input memory=%f sec, observation input memory=%f sec", stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i]))
        if(!is.nan(stateInputMemorySecs[i])) {
            stateInputs <- buildGoNogoAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=stateInputMemorySecs[i]*sRate)
        } else {
            stateInputs <- NA
        }
        if(!is.nan(obsInputMemorySecs[i])) {
            obsInputs <- buildGoNogoAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=obsInputMemorySecs[i]*sRate)
        } else {
            obsInputs <- NA
        }
        kem <- fit_FA_MARSS(observations=tSpikeRates, stateDim=stateDim[i], stateInputs=stateInputs, stateOffsetType=stateOffsetType, stateCovType=stateCovType, obsInputs=obsInputs, obsOffsetType=obsOffsetType, obsCovType=obsCovType, initialStateMeanType=initialStateMeanType, initialStateCovType=initialStateCovType, maxIter=maxIter, kfFunc=kfFunc)
        show(sprintf("number of latents=%d, state input memory=%f sec, observation input memory=%f sec: logLik=%f, AIC=%f, AICc=%f", stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i], kem$logLik, kem$AIC, kem$AICc))
        #
        resultsFilename <- sprintf(resultsFilenamePattern, stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i])
        results <- list(kem=kem, tSpikeRates=tSpikeRates, stateInputs=stateInputs, obsInputs=obsInputs, sRate=sRate, startTime=analysisStart)
        save(results, file=resultsFilename)
        return()
    }
}

processAll()
