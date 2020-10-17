
require(MASS)
require(ini)
require(optparse)
source("../projectSrc/utils/kalmanFilter/buildGoNogoVisualAndLaserInputs.R")
source("../commonSrc/stats/kalmanFilter/emEstimationKF_SS_withOffsetsAndInputs.R")
source("../commonSrc/stats/kalmanFilter/filterLDS_SS_withOffsetsAndInputs.R")
source("../commonSrc/stats/kalmanFilter/smoothLDS_SS.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondFA.R")
source("../commonSrc/stats/kalmanFilter/estimateKFInitialCondPPCA.R")

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

} else {
    # begin uncomment to debug
    stateDim <- 9
    stateInputMemorySecs <- 0.0
    obsInputMemorySecs <- 0.6
    initialCondMethod <- "FA"
    estConfigFilename <- "data/v1Shaft1_estimation_myEM.ini"
    modelsLogFilename <- "log/v1Shaft1_myEM.csv"
    # end uncomment to debug
}

    estConfig <- read.ini(estConfigFilename)

    stateCovType <- estConfig$covariance_type$states
    initialStateCovType <- estConfig$covariance_type$initialStates
    obsCovType <- estConfig$covariance_type$observations
    covsConstraints <- list(V0=initialStateCovType, Q=stateCovType, R=obsCovType)

    V0 <- eval(parse(text=estConfig$initial_values$V0))
    u0 <- eval(parse(text=estConfig$initial_values$u0))
    C0 <- eval(parse(text=estConfig$initial_values$C0))
    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    a0 <- eval(parse(text=estConfig$initial_values$a0))
    D0 <- eval(parse(text=estConfig$initial_values$D0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))

    M <- nrow(V0)

    if(tolower(estConfig$initial_values$m0)=="randomuniform") {
        m0Min <- as.double(estConfig$initial_value$m0Min)
        m0Max <- as.double(estConfig$initial_value$m0Max)
        m0 <- runif(n=M, min=m0Min, max=m0Max)
    } else {
        m0 <- eval(parse(text=estConfig$initial_values$m0))
    }

    region <- estConfig$control_variables$region
    shaft <- as.numeric(estConfig$control_variables$shaft)
    analysisStartTimeSecs <- as.double(estConfig$control_variables$analysisStartTimeSecs)
    analysisDurSecs <- as.double(estConfig$control_variables$analysisDurSecs)

    tol <- as.double(estConfig$EM$tol)
    maxIter <- as.numeric(estConfig$EM$maxIter)
    maxIter <- maxIter+1

    dataFilename <-  estConfig$filenames$dataFilename
    estMetaDataFilenamePattern <- estConfig$filenames$estMetaDataFilenamePattern
    estResFilenamePattern <- estConfig$filenames$estResFilenamePattern

    data <- get(load(dataFilename))
    sRate <- data$sRate
    analysisSamples <- (analysisStartTimeSecs*sRate)+(1:(analysisDurSecs*sRate))
    spikeCounts <- data[[sprintf("%sShaft%dSpikeCounts", region, shaft)]][,analysisSamples]
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

    metaData <- list()
    metaData[["estimation_config_info"]] <- list(estConfigFilename=estConfigFilename)
    metaData[["model_info"]] <- list(stateDim=stateDim, stateInputMem=stateInputMemorySecs, obsInputMem=obsInputMemorySecs)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    show(sprintf("Saved estimation meta data to: %s", estMetaDataFilename))

    startTime <- proc.time()[3]
    if(!is.nan(stateInputMemorySecs)) {
        stateInputs <- buildGoNogoVisualAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=as.integer(stateInputMemorySecs*sRate))
        dim(stateInputs) <- c(nrow(stateInputs), 1, ncol(stateInputs))
    } else {
        stateInputs <- NA
    }
    if(!is.nan(obsInputMemorySecs)) {
        obsInputs <- buildGoNogoVisualAndLaserInputs(goStim=goStim, nogoStim=nogoStim, laserStim=laserStim, inputMemory=as.integer(obsInputMemorySecs*sRate))
        dim(obsInputs) <- c(nrow(obsInputs), 1, ncol(obsInputs))
    } else {
        obsInputs <- NA
    }
    if(initialCondMethod=="FA") {
        dataForFA <- t(as.matrix(sqrtSpikeCounts))
        controlFA <- list(trace=TRUE, nstart=5)
        initialConds <- estimateKFInitialCondFA(z=dataForFA, nFactors=stateDim, control=controlFA)
        B0 <- initialConds$B
        Z0 <- initialConds$Z
        R0 <- diag(initialConds$RDiag)
        initialConds <- list(B0=B0, Z0=Z0, R0=R0, u0=u0, C0=C0, Q0=Q0, a0=a0, D0=D0, m0=m0, V0=V0)
    } else {
        if(initialCondMethod=="PPCA") {
            dataForPPCA <- t(as.matrix(sqrtSpikeCounts))
            initialConds <- estimateKFInitialCondPPCA(z=dataForPPCA, nFactors=stateDim)
            B0 <- initialConds$B
            Z0 <- initialConds$Z
            initialConds <- list(B=B0, Z=Z0, u=u0, C=C0, Q=Q0, a=a0, D=D0, R=R0, m0=m0, V0=V0)
        } else {
            stop(sprintf("Invalid initialCondMethod=%s", initialCondMethod))
        }
    }
    dsSSM <- emEstimationKF_SS_withOffsetsAndInputs(y=sqrtSpikeCounts, c=stateInputs, d=obsInputs, B0=initialConds$B, u0=initialConds$u, C0=initialConds$C, Q0=initialConds$Q, Z0=initialConds$Z, a0=initialConds$a, D0=initialConds$D, R0=initialConds$R, m0=initialConds$m0, V0=initialConds$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(m0=TRUE, V0=TRUE, B=TRUE, u=TRUE, C=TRUE, Q=TRUE, Z=TRUE, a=TRUE, D=TRUE, R=TRUE), covsConstraints=covsConstraints)
    elapsedTime <- proc.time()[3]-startTime
    elapsedTime <- unname(elapsedTime)
    AIC <- computeAIC(dsSSM=dsSSM)
    logMessage <- sprintf("%d, %d, %f, %f, %s, %f, %f, %f\n", estNumber, stateDim, stateInputMemorySecs, obsInputMemorySecs, initialCondMethod, dsSSM$logLik[length(dsSSM$logLik)], AIC, elapsedTime)
    show(logMessage)
    cat(logMessage, file=modelsLogFilename, append=TRUE)
    metaData[["estimation_summary"]] <- list(logLik=dsSSM$logLik[length(dsSSM$logLik)], AIC=AIC, elapsedTime=elapsedTime)
    write.ini(x=metaData, filepath=estMetaDataFilename)
    #
    estResFilename <- sprintf(estResFilenamePattern, estNumber)
    estRes <- list(dsSSM=dsSSM, AIC=AIC, initialConds=initialConds, stateDim=stateDim, stateInputMemorySecs=stateInputMemorySecs, obsInputMemorySecs=obsInputMemorySecs, sqrtSpikeCounts=sqrtSpikeCounts, stateInputs=stateInputs, obsInputs=obsInputs, sRate=sRate, startTime=analysisStartTimeSecs)
    save(estRes, file=estResFilename)
}

null <- processAll()
