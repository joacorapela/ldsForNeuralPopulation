
require(MASS)
require(MARSS)
require(ini)
require(optparse)
source("../projectSrc/utils/kalmanFilter/buildGoNogoVisualAndLaserInputs.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondPPCA.R")
source("../commonSrc/stats/kalmanFilter/fit_MARSS.R")
source("../commonSrc/stats/kalmanFilter/create_MARSS.R")

processAll <- function() {
DEBUG <- TRUE
if(!DEBUG) {
    option_list <- list( 
        make_option(c("-d", "--stateDim"), type="integer", default=3, help="State dimensionalty"),
        make_option(c("-s", "--stateInputMemorySecs"), type="double", default=0.6, help="State input memory (sec)"),
        make_option(c("-o", "--obsInputMemorySecs"), type="double", default=0.6, help="Observations input memory (sec)"),
        make_option(c("-i", "--initialCondMethod"), type="character", default="FA", help="Initial conditions method (FA: factor analysis; PPCA: probabilisitc PCA")
    )
    parser <- OptionParser(usage = "%prog [options] configFilename modelsLogFilename", option_list=option_list)
    parseRes <- parse_args(parser, positional_arguments=2)
    arguments <- parseRes$args
    options <- parseRes$options

    stateDim <- options$stateDim
    stateInputMemorySecs <- options$stateInputMemorySecs
    if(stateInputMemorySecs<0) {
        stateInputMemorySecs <- NaN
    }
    obsInputMemorySecs <- options$obsInputMemorySecs
    if(obsInputMemorySecs<0) {
        obsInputMemorySecs <- NaN
    }
    initialCondMethod <- options$initialCondMethod

    estConfigFilename <- arguments[[1]]
    modelsLogFilename <- arguments[[2]]

} else{
    # begin uncomment to debug
    stateDim <- 9
    stateInputMemorySecs <- 0.0
    obsInputMemorySecs <- 0.6
    initialCondMethod <- "FA"
    estConfigFilename <- "data/v1Shaft1_estimation.ini"
    modelsLogFilename <- "log/v1Shaf1lModels.csv"
    # end uncomment to debug
}

    estConfig <- read.ini(estConfigFilename)

    region <- estConfig$control_variables$region
    shafts <- eval(parse(text=estConfig$control_variables$shafts))
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

    data <- get(load(dataFilename))
    sRate <- data$sRate
    analysisSamples <- (analysisStartTimeSecs*sRate)+(1:(analysisDurSecs*sRate))
    spikeCounts <- c()
    for(shaft in shafts) {
        shaftSpikeCounts <- data[[sprintf("%sShaft%dSpikeCounts", region, shaft)]][,analysisSamples]
        nNeurons <- nrow(shaftSpikeCounts)
        rownames(shaftSpikeCounts) <- sprintf("shaft%02dNeuron%03d", rep(shaft, times=nNeurons), 1:nNeurons)
        spikeCounts <- rbind(spikeCounts, shaftSpikeCounts)
    }
    goStim <- data$goStim[analysisSamples]
    nogoStim <- data$nogoStim[analysisSamples]
    laserStim <- data$laserStim[analysisSamples]
    sqrtSpikeCounts <- sqrt(spikeCounts)

    show(sprintf("Processing number of latents=%d, state input memory=%f sec, observation input memory=%f sec", stateDim, stateInputMemorySecs, obsInputMemorySecs))

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
    metaData[["model_info"]] <- list(stateDim=stateDim, stateInputMem=stateInputMemorySecs, obsInputMem=obsInputMemorySecs)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    show(sprintf("Saved estimation meta data to: %s", estMetaDataFilename))

    startTime <- proc.time()[3]
    if(!is.nan(stateInputMemorySecs)) {
        stateInputs <- buildGoNogoVisualAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=as.integer(stateInputMemorySecs*sRate))
    } else {
        stateInputs <- NA
    }
    if(!is.nan(obsInputMemorySecs)) {
        obsInputs <- buildGoNogoVisualAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=as.integer(obsInputMemorySecs*sRate))
    } else {
        obsInputs <- NA
    }
    if(initialCondMethod=="FA") {
        dataForFA <- t(as.matrix(sqrtSpikeCounts))
        controlFA <- list(trace=TRUE, nstart=5)
        initialConds <- estimateKFInitialCondFA(z=dataForFA, nFactors=stateDim, control=controlFA)
        B0 <- matrix(as.vector(initialConds$B), ncol=1)
        Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
        R0 <- matrix(initialConds$RDiag, ncol=1)
        inits <- list(B=B0, Z=Z0, R=R0)
    } else {
        if(initialCondMethod=="PPCA") {
            dataForPPCA <- t(as.matrix(sqrtSpikeCounts))
            initialConds <- estimateKFInitialCondPPCA(z=dataForPPCA, nFactors=stateDim)
            B0 <- matrix(as.vector(initialConds$B), ncol=1)
            Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
            inits <- list(B=B0, Z=Z0)
        } else {
            stop(sprintf("Invalid initialCondMethod=%s", initialCondMethod))
        }
    }
    kem <- fit_MARSS(observations=sqrtSpikeCounts, inits=inits, stateDim=stateDim, stateInputs=stateInputs, stateOffsetType=stateOffsetType, stateCovType=stateCovType, obsInputs=obsInputs, obsOffsetType=obsOffsetType, obsCovType=obsCovType, initialStateMeanType=initialStateMeanType, initialStateCovType=initialStateCovType, maxIter=maxIter, kfFunc=kfFunc)
    elapsedTime <- proc.time()[3]-startTime
    kem <- MARSSaic(kem)
    logMessage <- sprintf("%d, %d, %f, %f, %s, %f, %f, %f, %f\n", estNumber, stateDim, stateInputMemorySecs, obsInputMemorySecs, initialCondMethod, kem$logLik, kem$AIC, kem$AICc, elapsedTime$elapsed)
    show(logMessage)
    cat(logMessage, file=modelsLogFilename, append=TRUE)
    metaData[["estimation_summary"]] <- list(logLik=kem$logLik, AIC=kem$AIC, AICc=kem$AICc, elapsedTime=elapsedTime)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    #
    estResFilename <- sprintf(estResFilenamePattern, estNumber)
    estRes <- list(kem=kem, stateDim=stateDim, stateInputMemorySecs=stateInputMemorySecs, obsInputMemorySecs=obsInputMemorySecs, sqrtSpikeCounts=sqrtSpikeCounts, stateInputs=stateInputs, obsInputs=obsInputs, sRate=sRate, startTime=analysisStartTimeSecs)
    save(estRes, file=estResFilename)
}

null <- processAll()
