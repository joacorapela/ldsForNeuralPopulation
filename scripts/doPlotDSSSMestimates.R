
require(MARSS)
require(ini)
require(optparse)
require(plotly)
source("../commonSrc/stats/utils/computePercentageExplainedVar.R")
source("../commonSrc/stats/kalmanFilter/filterLDS_SS_withOffsetsAndInputs.R")
source("../commonSrc/stats/kalmanFilter/smoothLDS_SS_withOffsetsAndInputs.R")
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
    DEBUG <- TRUE
    if(!DEBUG) {
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
    } else {
        estNumber <- 87430650
        estMetaDataFilenamePattern <- "results/%08d_estimation.ini"
        figFilenamePattern <- "figures/%08d_%s.%s"
    }
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

    dsSSM <- estRes$dsSSM
    obsInputs <- estRes$obsInputs
    stateInputs <- estRes$stateInputs
    sRate <- estRes$sRate
    stateInputMemorySamples <- stateInputMemorySecs*sRate
    obsInputMemorySamples <- obsInputMemorySecs*sRate

    kfRes <- filterLDS_SS_withOffsetsAndInputs(y=estRes$sqrtSpikeCounts, B=dsSSM$B, u=dsSSM$u, C=dsSSM$C, c=estRes$stateInputs, Q=dsSSM$Q, m0=dsSSM$m0, V0=dsSSM$V0, Z=dsSSM$Z, a=dsSSM$a, D=dsSSM$D, d=estRes$obsInputs, R=dsSSM$R)
    ksRes <- smoothLDS_SS_withOffsetsAndInputs(B=dsSSM$B, u=dsSSM$u, C=dsSSM$C, c=dsSSM$c, Q=dsSSM$Q, xnn=kfRes$xnn, Vnn=kfRes$Vnn, xnn1=kfRes$xnn1, Vnn1=kfRes$Vnn1, m0=dsSSM$m0, V0=dsSSM$V0, initStateAt=0)

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

if(FALSE) {
    show("Plotting logLik")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "html")
    fig <- getPlotLogLik(logLik=dsSSM$logLik)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting B")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "B", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "B", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=dsSSM$B)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting stateSpectrum")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "html")
    fig <- getPlotStateSpectrum(B=dsSSM$B)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting U")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "U", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "U", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$u)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(stateInputMemorySecs)) {
        show("Plotting C")
        x <- (0:stateInputMemorySamples)/sRate
        CblockSize <- 1+stateInputMemorySamples
        Coffset <- 0
        #
        CVisualStimGo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CVisualStimNogo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimGo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimNogo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting Q")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Q", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Q", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=diag(dsSSM$Q), xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting m0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "m0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "m0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$m0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting V0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "V0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "V0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=diag(dsSSM$V0), xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting Z")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Z", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Z", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=dsSSM$Z)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting a")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "a", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "a", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$a, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(obsInputMemorySecs)) {
        show("Plotting D")
        x <- (0:obsInputMemorySamples)/sRate
        DblockSize <- 1+obsInputMemorySamples
        Doffset <- 0
        #
        DVisualStimGo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DVisualStimGo, estimatedLegendLabelPattern="neuron %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        DVisualStimNogo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DVisualStimNogo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        DLaserStimGo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DLaserStimGo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset+DblockSize
        #
        DLaserStimNogo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DLaserStimNogo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting R")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "R", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "R", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$R, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

}

    predStats <- computeOneStepAheadObsPredStats(xtt1=kfRes$xnn1[,1,], Vtt1=kfRes$Vnn1, Z=dsSSM$Z, a=as.vector(dsSSM$a), D=dsSSM$D, R=dsSSM$R, obsInputs=estRes$obsInputs[,1,])

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
    fig <- getPlotSmoothedStates(time=time, xtT=ksRes$xnN, VtT=ksRes$VnN, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)


    for(i in 1:nrow(ksRes$xnN)) {
        show(sprintf("Plotting smoothedState %d", i))
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=ksRes$xnN, VtT=ksRes$VnN, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff, statesToPlot=c(i))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    browser()
}

processAll()
