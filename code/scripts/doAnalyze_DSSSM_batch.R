
require(MASS)
require(ini)
require(optparse)
source("../projectSrc/utils/kalmanFilter/buildGoNogoVisualAndLaserInputs.R")
source("../projectSrc/utils/kalmanFilter/getInitialConds.R")
source("../commonSrc/stats/kalmanFilter/emEstimationKF_SS_withOffsetsAndInputs.R")
source("../commonSrc/stats/kalmanFilter/filterLDS_SS_withOffsetsAndInputs.R")
source("../commonSrc/stats/kalmanFilter/smoothLDS_SS.R")
source("../commonSrc/stats/kalmanFilter/computeAIC.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondPPCA.R")

modelInLog <- function(modelsLogFilename, analysisStartTimeSecs, trainDurSecs, validationDurSecs, stateDim, stateInputMemorySecs, obsInputMemorySecs, initialCondMethod) {
    # logMessage <- sprintf("%d, %f, %f, %f, %d, %f, %f, %s, %f, %f, %f, %f\n", estNumber, analysisStartTimeSecs, trainDurSecs, validationDurSecs, stateDim, stateInputMemorySecs, obsInputMemorySecs, initialCondMethod, dsSSM$logLik[length(dsSSM$logLik)], AIC, cvLogLike, elapsedTime)
    if(!files.exists(modelsLogFilename)) {
           return FALSE
    }

    lines <- readLines(con=modelsLogFilename)
    tokens <- strsplit(lines, ",")
    found <- FALSE
    i <- 1
    while(i<=length(tokens) && !found) {
        tAnalysisStartTimeSecs <- as.numeric(tokens[[i]][2])
        tTrainDurSecs <- as.numeric(tokens[[i]][3])
        tValidationDurSecs <- as.numeric(tokens[[i]][4])
        tStateDim <- as.numeric(tokens[[i]][5])
        tStateInputMemSecs <- as.numeric(tokens[[i]][6])
        tObsInpuMemSecs <- as.numeric(tokens[[i]][7])
        tInitialCondMethod <- gsub(" ", "", tokens[[i]][8], fixed=TRUE)
        if(analysisStartTimeSecs==tAnalysisStartTimeSecs && trainDurSecs==tTrainDurSecs && validationDurSecs==tValidationDurSecs && stateDim==tStateDim && stateInputMemorySecs==tStateInputMemSecs && obsInputMemorySecs==tObsInpuMemSecs && initialCondMethod==tInitialCondMethod) {
            found <- TRUE
        }
        i <- i+1
    }
    return(found)
}

processAll <- function() {
# DEBUG <- TRUE
DEBUG <- FALSE
if(!DEBUG) {
    option_list <- list( 
        make_option(c("-a", "--analysisStartTimeSecs"), type="double", default=0, help="Analysis start time (sec)"),
        make_option(c("-t", "--trainDurSecs"), type="double", default=180, help="Duration of the data train segment (sec)"),
        make_option(c("-v", "--validationDurSecs"), type="double", default=60, help="Duration of the data validation segment (sec)"),
        make_option(c("-d", "--stateDim"), type="integer", default=3, help="State dimensionalty"),
        make_option(c("-m", "--stateInputMemorySecs"), type="double", default=0.6, help="State input memory (sec)"),
        make_option(c("-o", "--obsInputMemorySecs"), type="double", default=0.6, help="Observations input memory (sec)"),
        make_option(c("-i", "--initialCondMethod"), type="character", default="FA", help="Initial conditions method (FA: factor analysis; PPCA: probabilisitc PCA"),
        make_option(c("-n", "--nStartFA"), type="integer", default=5, help="Number of start values for factor analysis"),
        make_option(c("-c", "--checkModelsLogFilename"), action="store_true", default=FALSE, help="Check models log filename and do not estimate a model if it has been previously estimated")
    )
    parser <- OptionParser(usage = "%prog [options] binConfigFilename estConfigFilename modelsLogFilename", option_list=option_list)
    parseRes <- parse_args(parser, positional_arguments=3)
    arguments <- parseRes$args
    options <- parseRes$options

    analysisStartTimeSecs <- options$analysisStartTimeSecs
    trainDurSecs <- options$trainDurSecs
    validationDurSecs <- options$validationDurSecs

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
    nStartFA <- options$nStartFA
    checkModelsLogFilename <- options$checkModelsLogFilename

    binConfigFilename <- arguments[[1]]
    estConfigFilename <- arguments[[2]]
    modelsLogFilename <- arguments[[3]]

} else {
    # 46390734, 180.000000, 180.000000, 60.000000, 5, 0.000000, 0.400000, PPCA, -10923.367918, 23110.735836, -3692.318680, 141.450000
    analysisStartTimeSecs <- 180
    trainDurSecs <- 180
    validationDurSecs <- 60
    stateDim <- 50
    stateInputMemorySecs <- 0.0
    obsInputMemorySecs <- 0.4
    initialCondMethod <- "PPCA"
    binConfigFilename <- "../../data/VL61/binLDStimeSeries.ini"
    estConfigFilename <- "../../data/VL61/v1Shaft1_estimation_DSSSM.ini"
    modelsLogFilename <- "../../log/VL61/v1Shaft1Models_DSSSM.csv"
    # estConfigFilename <- "data/v1Shaft1_estimation_DSSSM_start2580_dur600.ini"
    # modelsLogFilename <- "log/v1Shaft1Models_DSSSM_start2580_dur600.csv"
    nStartFA <- 5
    checkModelsLogFilename <- TRUE
}
    if(modelInLog(modelsLogFilename=modelsLogFilename, analysisStartTimeSecs=analysisStartTimeSecs, trainDurSecs=trainDurSecs, validationDurSecs=validationDurSecs, stateDim=stateDim, stateInputMemorySecs=stateInputMemorySecs, obsInputMemorySecs=obsInputMemorySecs, initialCondMethod=initialCondMethod)) {
           return()
    }

    binConfig <- read.ini(binConfigFilename)
    dataFilename <-  binConfig$filenames$saveFilename

    estConfig <- read.ini(estConfigFilename)

    stateCovType <- estConfig$covariance_type$states
    initialStateCovType <- estConfig$covariance_type$initialStates
    obsCovType <- estConfig$covariance_type$observations
    covsConstraints <- list(V0=initialStateCovType, Q=stateCovType, R=obsCovType)

    region <- estConfig$control_variables$region
    shafts <- eval(parse(text=estConfig$control_variables$shafts))

    tol <- as.double(estConfig$EM$tol)
    maxIter <- as.numeric(estConfig$EM$maxIter)
    maxIter <- maxIter+1

    estMetaDataFilenamePattern <- estConfig$filenames$estMetaDataFilenamePattern
    estResFilenamePattern <- estConfig$filenames$estResFilenamePattern

    data <- get(load(dataFilename))
    sRate <- data$sRate
    trainSamples <- (analysisStartTimeSecs*sRate)+(1:(trainDurSecs*sRate))
    validationSamples <- max(trainSamples)+(1:(validationDurSecs*sRate))

    spikeCounts <- c()
    for(shaft in shafts) {
        shaftSpikeCounts <- data[[sprintf("%sShaft%dSpikeCounts", region, shaft)]]
        nNeurons <- nrow(shaftSpikeCounts)
        rownames(shaftSpikeCounts) <- sprintf("shaft%02dNeuron%03d", rep(shaft, times=nNeurons), 1:nNeurons)
        spikeCounts <- rbind(spikeCounts, shaftSpikeCounts)
    }
    trainSpikeCounts <- spikeCounts[,trainSamples]
    trainSqrtSpikeCounts <- sqrt(trainSpikeCounts)
    validationSpikeCounts <- spikeCounts[,validationSamples]
    validationSqrtSpikeCounts <- sqrt(validationSpikeCounts)

    trainGoStim <- data$goStim[trainSamples]
    validationGoStim <- data$goStim[validationSamples]

    trainNogoStim <- data$nogoStim[trainSamples]
    validationNogoStim <- data$nogoStim[validationSamples]

    trainLaserStim <- data$laserStim[trainSamples]
    validationLaserStim <- data$laserStim[validationSamples]

    obsDim <- nrow(trainSqrtSpikeCounts)
    initialConds <- getInitialConds(initialValues=estConfig$initial_values, stateDim=stateDim, obsDim=obsDim, stateInputMemorySamples=stateInputMemorySecs*sRate, obsInputMemorySamples=obsInputMemorySecs*sRate)

    show(sprintf("Processing number of latents=%d, state input memory=%f sec, observation input memory=%f sec", stateDim, stateInputMemorySecs, obsInputMemorySecs))

    exit <- FALSE
    while(!exit) {
        estNumber <- sample(1e8, 1)
        estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
        if(!file.exists(estMetaDataFilename)) {
            exit <- TRUE
        }
    }

    metaData <- list()
    metaData[["estimation_config_info"]] <- list(estConfigFilename=estConfigFilename)
    metaData[["model_info"]] <- list(stateDim=stateDim, stateInputMem=stateInputMemorySecs, obsInputMem=obsInputMemorySecs)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    show(sprintf("Saved estimation meta data to: %s", estMetaDataFilename))

    if(!is.nan(stateInputMemorySecs)) {
        trainStateInputs <- buildGoNogoVisualAndLaserInputs(goStim=trainGoStim, nogoStim=trainNogoStim, laserStim=trainLaserStim, inputMemory=as.integer(stateInputMemorySecs*sRate))
        dim(trainStateInputs) <- c(nrow(trainStateInputs), 1, ncol(trainStateInputs))
        validationStateInputs <- buildGoNogoVisualAndLaserInputs(goStim=validationGoStim, nogoStim=validationNogoStim, laserStim=validationLaserStim, inputMemory=as.integer(stateInputMemorySecs*sRate))
        dim(validationStateInputs) <- c(nrow(validationStateInputs), 1, ncol(validationStateInputs))
    } else {
        trainStateInputs <- NA
        validationStateInputs <- NA
    }
    if(!is.nan(obsInputMemorySecs)) {
        trainObsInputs <- buildGoNogoVisualAndLaserInputs(goStim=trainGoStim, nogoStim=trainNogoStim, laserStim=trainLaserStim, inputMemory=as.integer(obsInputMemorySecs*sRate))
        dim(trainObsInputs) <- c(nrow(trainObsInputs), 1, ncol(trainObsInputs))
        validationObsInputs <- buildGoNogoVisualAndLaserInputs(goStim=validationGoStim, nogoStim=validationNogoStim, laserStim=validationLaserStim, inputMemory=as.integer(obsInputMemorySecs*sRate))
        dim(validationObsInputs) <- c(nrow(validationObsInputs), 1, ncol(validationObsInputs))
    } else {
        trainObsInputs <- NA
        validationObsInputs <- NA 
    }

    dataForEstInitialCond <- t(as.matrix(trainSqrtSpikeCounts))
    startTime <- proc.time()[3]
    if(initialCondMethod=="FA") {
        controlFA <- list(trace=TRUE, nstart=nStartFA)
        estRes <- estimateKFInitialCondFA(z=dataForEstInitialCond, nFactors=stateDim, control=controlFA)
        initialConds <- c(list(B=estRes$B, Z=estRes$Z, R=diag(estRes$RDiag)), initialConds)
    } else {
        if(initialCondMethod=="PPCA") {
            estRes <- estimateKFInitialCondPPCA(z=dataForEstInitialCond, nFactors=stateDim)
            initialConds <- c(list(B=estRes$B, Z=estRes$Z), initialConds)
        } else {
            stop(sprintf("Invalid initialCondMethod=%s", initialCondMethod))
        }
    }
    dsSSM <- emEstimationKF_SS_withOffsetsAndInputs(y=trainSqrtSpikeCounts, c=trainStateInputs, d=trainObsInputs, B0=initialConds$B, u0=initialConds$u, C0=initialConds$C, Q0=initialConds$Q, Z0=initialConds$Z, a0=initialConds$a, D0=initialConds$D, R0=initialConds$R, m0=initialConds$m0, V0=initialConds$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(m0=TRUE, V0=TRUE, B=TRUE, u=TRUE, C=TRUE, Q=TRUE, Z=TRUE, a=TRUE, D=TRUE, R=TRUE), covsConstraints=covsConstraints)
    elapsedTime <- proc.time()[3]-startTime
    elapsedTime <- unname(elapsedTime)

    AIC <- computeAIC(dsSSM=dsSSM)
    cvLogLike <- filterLDS_SS_withOffsetsAndInputs(y=validationSqrtSpikeCounts, c=validationStateInputs, d=validationObsInputs, B=dsSSM$B, u=dsSSM$u, C=dsSSM$C, Q=dsSSM$Q, m0=dsSSM$xNN, V0=dsSSM$VNN, Z=dsSSM$Z, a=dsSSM$a, D=dsSSM$D, R=dsSSM$R)$logLike
    logMessage <- sprintf("%d, %f, %f, %f, %d, %f, %f, %s, %f, %f, %f, %f\n", estNumber, analysisStartTimeSecs, trainDurSecs, validationDurSecs, stateDim, stateInputMemorySecs, obsInputMemorySecs, initialCondMethod, dsSSM$logLik[length(dsSSM$logLik)], AIC, cvLogLike, elapsedTime)
    show(logMessage)
    cat(logMessage, file=modelsLogFilename, append=TRUE)
    metaData[["estimation_summary"]] <- list(logLik=dsSSM$logLik[length(dsSSM$logLik)], AIC=AIC, cvLogLike=cvLogLike, elapsedTime=elapsedTime)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    #
    estResFilename <- sprintf(estResFilenamePattern, estNumber)
    estRes <- list(dsSSM=dsSSM, AIC=AIC, cvLogLike=cvLogLike, initialConds=initialConds, stateDim=stateDim, stateInputMemorySecs=stateInputMemorySecs, obsInputMemorySecs=obsInputMemorySecs, trainSqrtSpikeCounts=trainSqrtSpikeCounts, stateInputs=trainStateInputs, obsInputs=trainObsInputs, sRate=sRate, startTime=analysisStartTimeSecs)
    save(estRes, file=estResFilename)
}

null <- processAll()
