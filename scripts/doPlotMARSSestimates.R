
require(MARSS)
require(ini)
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
    estNumber <- 33778847 # v1Shaft1estMetaData
    estMetaDataFilenamePattern <- "results/%08d_estimation.ini"
    estConfigFilenamePattern <- "data/%08d_estimation.ini"
    figFilenamePattern <- "figures/%08d_%s.%s"

    estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
    estMetaData <- read.ini(filepath=estMetaDataFilename)

    estConfigFilename <- estMetaData$estimation_config_info$estConfigFilename
    estConfig <- read.ini(estConfigFilename)
    # obsInputMemorySecs <- as.double(estConfig$inputs$obsInputMemorySecs)
    # stateInputMemorySecs <- as.double(estConfig$inputs$stateInputMemorySecs)
    # stateDim <- eval(parse(text=estConfig$control_variables$stateDim))
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
#     kem0 <- kem
#     kem0$par <- kem0$start
#     kfRes0 <- MARSSkf(kem0)

    sqrtSpikeCounts <- estRes$sqrtSpikeCounts
    time <- estRes$startTime+(1:ncol(estRes$sqrtSpikeCounts))/estRes$sRate
    dimObs <- nrow(sqrtSpikeCounts)
    nObs <- ncol(sqrtSpikeCounts)

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
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$B, estimated=matrix(coef(kem)$B, nrow=stateDim))
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

    show("Plotting C")
    if(!is.nan(stateInputMemorySecs)) {
        x <- (0:stateInputMemorySamples)/sRate
        CblockSize <- stateDim*(1+stateInputMemorySamples)
        Coffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Cgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Cgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim), estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Cnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Cnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim), estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Claser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Claser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim), estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CgoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CgoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim), estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CnogoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CnogoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim), estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
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
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$Z, estimated=matrix(coef(kem)$Z, nrow=dimObs))
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

    show("Plotting D")
    if(!is.nan(obsInputMemorySecs)) {
        x <- (0:obsInputMemorySamples)/sRate
        DblockSize <- dimObs*(1+obsInputMemorySamples)
        Doffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dlaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dlaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DgoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DgoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs), estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DnogoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DnogoLaser", "html")
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
    fig <- getPlotOneStepAheadForecasts(time=time, obs=sqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, obsInputMemorySamples=obsInputMemorySamples, inputs=obsInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting oneStepAheadForecastsNeuron")
    for(i in 1:nrow(predStats$ytt1)) {
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "html")
        fig <- getPlotOneStepAheadForecasts(time=time, obs=sqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, obsToPlot=c(i), obsInputMemorySamples=obsInputMemorySamples, inputs=obsInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting smoothedStates")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "html")
    fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, stateInputMemorySamples=stateInputMemorySamples, inputs=stateInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    for(i in 1:nrow(kfRes$xtT)) {
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smoothedState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, statesToPlot=c(i), stateInputMemorySamples=stateInputMemorySamples, inputs=stateInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    browser()
}

processAll()
