
require(MARSS)
require(ini)
require(optparse)
require(plotly)
source("../commonSrc/stats/utils/computePercentageExplainedVar.R")
source("../commonSrc/stats/kalmanFilter/computeOneStepAheadObsPredStats.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedMatrices.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedVectors.R")
source("../commonSrc/plot/kalmanFilter/getPlotStateSpectrum.R")
source("../commonSrc/plot/kalmanFilter/getStimRecs.R")
source("../commonSrc/plot/kalmanFilter/getPlotOneStepAheadForecasts.R")
source("../commonSrc/plot/kalmanFilter/getPlotSmoothedStates.R")
source("../commonSrc/plot/kalmanFilter/getPlotPercentageExplainedVar.R")
source("../commonSrc/plot/kalmanFilter/getPlotLogLik.R")

processAll <- function() {
    option_list <- list( 
        make_option(c("-m", "--estMetaDataFilenamePattern"), type="character", default="results/%08d_estimation.ini", help="Estimation metadata filename pattern"),
        make_option(c("-f", "--figFilenamePattern"), type="character", default="figures/%08d_%s.%s", help="Figure filename pattern")
    )
    parser <- OptionParser(usage = "%prog [options] estNumber", option_list=option_list)
    parseRes <- parse_args(parser, positional_arguments=1)
    options <- parseRes$options
    arguments <- parseRes$args

    estNumber <- as.numeric(arguments[1])
    estMetaDataFilenamePattern <- options$estMetaDataFilenamePattern
    figFilenamePattern <- options$figFilenamePattern

    estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
    estMetaData <- read.ini(filepath=estMetaDataFilename)
    estConfigFilename <- estMetaData$estimation_config_info$estConfigFilename
    estConfig <- read.ini(estConfigFilename)
    dataFilename <- estConfig$filenames$dataFilename
    estResFilenamePattern <- estConfig$filenames$estResFilenamePattern

    estResFilename <- sprintf(estResFilenamePattern, estNumber)
    estRes <- get(load(estResFilename))
    stateDim <- estRes$stateDim
    obsInputMemorySecs <- estRes$obsInputMemorySecs
    stateInputMemorySecs <- estRes$stateInputMemorySecs

    kem <- estRes$kem
    obsInputs <- estRes$obsInputs
    stateInputs <- estRes$stateInputs
    initialConds <- list(B=matrix(kem$start$B, ncol=stateDim), Z=matrix(kem$start$Z, ncol=stateDim), RDiag=as.vector(kem$start$R))
    sRate <- estRes$sRate
    stateInputMemorySamples <- stateInputMemorySecs*sRate
    obsInputMemorySamples <- obsInputMemorySecs*sRate

    kfRes <- MARSSkf(kem)

    sqrtSpikeCounts <- estRes$sqrtSpikeCounts
    time <- estRes$startTime+(1:ncol(estRes$sqrtSpikeCounts))/estRes$sRate
    dimObs <- nrow(sqrtSpikeCounts)
    nObs <- ncol(sqrtSpikeCounts)

    data <- get(load(dataFilename))
    minTime <- min(time)
    maxTime <- max(time)

    goStimOn <- data$goStimOnSecs
    goStimOff <- data$goStimOffSecs
    toKeepIndices <- which(minTime<=goStimOn & goStimOff<=maxTime)
    goStimOn <- goStimOn[toKeepIndices]
    goStimOff <- goStimOff[toKeepIndices]

    nogoStimOn <- data$nogoStimOnSecs
    nogoStimOff <- data$nogoStimOffSecs
    toKeepIndices <- which(minTime<=nogoStimOn & nogoStimOff<=maxTime)
    nogoStimOn <- nogoStimOn[toKeepIndices]
    nogoStimOff <- nogoStimOff[toKeepIndices]

    laserStimOn <- data$laserStimOnSecs
    laserStimOff <- data$laserStimOffSecs
    toKeepIndices <- which(minTime<=laserStimOn & laserStimOff<=maxTime)
    laserStimOn <- laserStimOn[toKeepIndices]
    laserStimOff <- laserStimOff[toKeepIndices]

# if(FALSE) {
    show("Plotting logLik")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "html")
    fig <- getPlotLogLik(logLik=estRes$kem$iter.record$logLik)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting B")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "B", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "B", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$B, nrow=stateDim))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting stateSpectrum")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "html")
    fig <- getPlotStateSpectrum(B=matrix(coef(kem)$B, nrow=stateDim))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting U")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "U", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "U", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$U, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(stateInputMemorySecs)) {
        show("Plotting C")
        x <- (0:stateInputMemorySamples)/sRate
        CblockSize <- stateDim*(1+stateInputMemorySamples)
        Coffset <- 0
        #
        CVisualStimGo <- matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CVisualStimNogo <- matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimGo <- matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimNogo <- matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
    }

    show("Plotting Q")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Q", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Q", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting m0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "m0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "m0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting V0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "V0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "V0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$V0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting Z")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Z", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Z", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$Z, nrow=dimObs))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting a")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "a", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "a", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$A, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(obsInputMemorySecs)) {
        show("Plotting D")
        x <- (0:obsInputMemorySamples)/sRate
        DblockSize <- dimObs*(1+obsInputMemorySamples)
        Doffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
    }

    show("Plotting R")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "R", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "R", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$R, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

# }

    if(kem$call$model$R=="diagonal and unequal") {
        R <- diag(coef(kem)$R)
    } else {
        error("Functionality not yet implemented for R!=<diagonal and unequal>")
    }
    predStats <- computeOneStepAheadObsPredStats(xtt1=kfRes$xtt1, Vtt1=kfRes$Vtt1, Z=matrix(coef(kem)$Z, nrow=dimObs), a=as.numeric(coef(kem)$A), D=matrix(coef(kem)$D, nrow=dimObs), R=R, obsInputs=obsInputs)

    show("Plotting percExpVar")
    percExpVar <- computePercentageExplainedVar(observations=sqrtSpikeCounts, predictions=predStats$ytt1)
    pngFilename <- sprintf(figFilenamePattern, estNumber, "percExpVar", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "percExpVar", "html")
    fig <- getPlotPercentageExplainedVar(percExpVar=percExpVar)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting oneStepAheadForecasts")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "oneStepAheadForecasts", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "oneStepAheadForecasts", "html")
    fig <- getPlotOneStepAheadForecasts(time=time, obs=sqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    for(i in 1:nrow(predStats$ytt1)) {
        show(sprintf("Plotting oneStepAheadForecast for neuron %d", i))
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "html")
        fig <- getPlotOneStepAheadForecasts(time=time, obs=sqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff, obsToPlot=c(i))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting smoothedStates")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "html")
    fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)


    for(i in 1:nrow(kfRes$xtT)) {
        show(sprintf("Plotting smoothedState %d", i))
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff, statesToPlot=c(i))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    browser()
}

processAll()
