
require(MASS)
require(MARSS)
require(ini)
source("../projectSrc/utils/kalmanFilter/buildGoNogoLaserAndInteractionsInputs.R")
source("../projectSrc/utils/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/fit_FA_MARSS.R")
source("../commonSrc/stats/kalmanFilter/create_FA_MARSS.R")

processAll <- function() {
    estConfigNumber <- 1
    estConfigFilenamePattern <- "data/%08d_estimation.ini"

    estConfigFilename <- sprintf(estConfigFilenamePattern, estConfigNumber)
    estConfig <- read.ini(estConfigFilename)

    obsInputMemorySecs <- eval(parse(text=estConfig$inputs$obsInputMemorySecs))
    stateInputMemorySecs <- eval(parse(text=estConfig$inputs$stateInputMemorySecs))

    region <- estConfig$control_variables$region
    shaft <- as.numeric(estConfig$control_variables$shaft)
    stateDim <- eval(parse(text=estConfig$control_variables$stateDim))
    analysisStartTimeSecs <- as.double(estConfig$control_variables$analysisStartTimeSecs)
    analysisDurSecs <- as.double(estConfig$control_variables$analysisDurSecs)

    maxIter <- as.numeric(estConfig$EM$maxIter)

    kfFunc <- estConfig$MARSS$kfFunc
    stateOffsetType <- estConfig$MARSS$stateOffsetType
    stateCovType <- estConfig$MARSS$stateCovType
    obsOffsetType <- estConfig$MARSS$obsOffsetType
    obsCovType <- estConfig$MARSS$obsCovType
    initialStateMeanType <- estConfig$MARSS$initialStateMeanType
    initialStateCovType <- estConfig$MARSS$initialStateCovType

    dataFilename <-  estConfig$filenames$dataFilename
    estMetaDataFilenamePattern <- estConfig$filenames$estMetaDataFilenamePattern
    estResFilenamePattern <- estConfig$filenames$estResFilenamePattern

    stopifnot(length(stateDim)==length(stateInputMemorySecs) & length(stateInputMemorySecs)==length(obsInputMemorySecs))
    #
    data <- get(load(dataFilename))
    sRate <- data$sRate
    binSizeSecs <- (data$breaks[2]-data$breaks[1])/sRate
    analysisSamples <- (analysisStartTimeSecs*sRate)+(1:(analysisDurSecs*sRate))
    spikeCounts <- data[[sprintf("%sShaft%dSpikeCounts", region, shaft)]][,analysisSamples]
    goStim <- data$goStim[analysisSamples]
    nogoStim <- data$nogoStim[analysisSamples]
    laserStim <- data$laserStim[analysisSamples]
    spikeRates <- spikeCounts/binSizeSecs
    tSpikeRates <- sqrt(spikeRates)

    for(i in 1:length(stateDim)) {
        show(sprintf("Processing number of latents=%d, state input memory=%f sec, observation input memory=%f sec", stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i]))

        exit <- FALSE
        while(!exit) {
            estNumber <- sample(1e8, 1)
            estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
            if(!file.exists(estMetaDataFilename)) {
                exit <- TRUE
            }
        }

        estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
        metaData <- list()
        metaData[["estimation_config_info"]] <- list(estConfigNumber=estConfigNumber)
        metaData[["model_info"]] <- list(stateDim=stateDim[i], stateInputMem=stateInputMemorySecs[i], obsInputMem=obsInputMemorySecs[i])
        write.ini(x=metaData, filepath=estMetaDataFilename)
        show(sprintf("Saved estimation meta data to: %s", estMetaDataFilename))

        if(!is.nan(stateInputMemorySecs[i])) {
            stateInputs <- buildGoNogoLaserAndInteractionsInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=stateInputMemorySecs[i]*sRate)
        } else {
            stateInputs <- NA
        }
        if(!is.nan(obsInputMemorySecs[i])) {
            obsInputs <- buildGoNogoLaserAndInteractionsInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=obsInputMemorySecs[i]*sRate)
        } else {
            obsInputs <- NA
        }
        kem <- fit_FA_MARSS(observations=tSpikeRates, stateDim=stateDim[i], stateInputs=stateInputs, stateOffsetType=stateOffsetType, stateCovType=stateCovType, obsInputs=obsInputs, obsOffsetType=obsOffsetType, obsCovType=obsCovType, initialStateMeanType=initialStateMeanType, initialStateCovType=initialStateCovType, maxIter=maxIter, kfFunc=kfFunc)
        kem <- MARSSaic(kem)
        show(sprintf("number of latents=%d, state input memory=%f sec, observation input memory=%f sec: logLik=%f, AIC=%f, AICc=%f", stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i], kem$logLik, kem$AIC, kem$AICc))
        #
        estResFilename <- sprintf(estResFilenamePattern, estNumber)
        estRes <- list(kem=kem, stateDim=stateDim[i], stateInputMem=stateInputMemorySecs[i], obsInputMem=obsInputMemorySecs[i], tSpikeRates=tSpikeRates, stateInputs=stateInputs, obsInputs=obsInputs, sRate=sRate, startTime=analysisStartTimeSecs)
        save(estRes, file=estResFilename)
    }
}

processAll()
