
require(MASS)
require(MARSS)
require(ini)
require(optparse)
source("../projectSrc/utils/kalmanFilter/buildGoNogoLaserAndInteractionsInputs.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondPPCA.R")
source("../commonSrc/stats/kalmanFilter/fit_MARSS.R")
source("../commonSrc/stats/kalmanFilter/create_MARSS.R")

processAll <- function() {
    option_list <- list( 
        make_option(c("-d", "--stateDim"), type="integer", default=3, help="State dimensionalty"),
        make_option(c("-s", "--stateInputMemorySecs"), type="double", default=0.6, help="State input memory (sec)"),
        make_option(c("-o", "--obsInputMemorySecs"), type="double", default=0.6, help="Observations input memory (sec)"),
        make_option(c("-i", "--initialCondMethod"), type="character", default="FA", help="Initial conditions method (FA: factor analysis; PPCA: probabilisitc PCA")
    )
    parser <- OptionParser(usage = "%prog [options] configFilename", option_list=option_list)
    parseRes <- parse_args(parser, positional_arguments=1)
    arguments <- parseRes$args
    options <- parseRes$options

    stateDim <- options$stateDim
    stateInputMemorySecs <- options$stateInputMemorySecs
    obsInputMemorySecs <- options$obsInputMemorySecs
    initialCondMethod <- options$initialCondMethod

    estConfigFilename <- arguments[[1]]
    estConfig <- read.ini(estConfigFilename)

    region <- estConfig$control_variables$region
    shaft <- as.numeric(estConfig$control_variables$shaft)
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
    modelsLogFilename <- estConfig$filenames$modelsLogFilename

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
        metaData[["estimation_config_info"]] <- list(estConfigFilename=estConfigFilename)
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
        if(initialCondMethod=="FA") {
            dataForFA <- t(as.matrix(tSpikeRates))
            initialConds <- estimateKFInitialCondFA(z=dataForFA, nFactors=stateDim, control=controlFA)
            B0 <- matrix(as.vector(initialConds$B), ncol=1)
            Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
            R0 <- matrix(initialConds$RDiag, ncol=1)
            inits <- list(B=B0, Z=Z0, R=R0)
        } else {
            if(initialCondMethod=="PPCA") {
                dataForPPCA <- t(as.matrix(tSpikeRates))
                initialConds <- estimateKFInitialCondPPCA(z=dataForPPCA, nFactors=stateDim)
                B0 <- matrix(as.vector(initialConds$B), ncol=1)
                Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
                inits <- list(B=B0, Z=Z0)
            } else {
                stop(sprintf("Invalid initialCondMethod=%s", initialCondMethod))
            }
        }
        kem <- fit_MARSS(observations=tSpikeRates, inits=inits, stateDim=stateDim[i], stateInputs=stateInputs, stateOffsetType=stateOffsetType, stateCovType=stateCovType, obsInputs=obsInputs, obsOffsetType=obsOffsetType, obsCovType=obsCovType, initialStateMeanType=initialStateMeanType, initialStateCovType=initialStateCovType, maxIter=maxIter, kfFunc=kfFunc)
        }
        kem <- MARSSaic(kem)
        logMessage <- sprintf("%d: number of latents=%d, state input memory=%f sec, observation input memory=%f sec, init=%s: logLik=%f, AIC=%f, AICc=%f", estNumber, stateDim[i], stateInputMemorySecs[i], obsInputMemorySecs[i], initialCondMethod, kem$logLik, kem$AIC, kem$AICc)
        writeLines(text=logMessage, con=modelsLogFilename)
        show(logMessage)
        metaData[["estimation_summary"]] <- list(logLik=kem$logLik, AIC=kem$AIC, AICc=kem$AICc)
        write.ini(x=metaData, filepath=estMetaDataFilename)

        #
        estResFilename <- sprintf(estResFilenamePattern, estNumber)
        estRes <- list(kem=kem, stateDim=stateDim[i], stateInputMem=stateInputMemorySecs[i], obsInputMem=obsInputMemorySecs[i], tSpikeRates=tSpikeRates, stateInputs=stateInputs, obsInputs=obsInputs, sRate=sRate, startTime=analysisStartTimeSecs)
        save(estRes, file=estResFilename)
    }
}

processAll()
