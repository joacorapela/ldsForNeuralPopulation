
require(MARSS)
require(plotly)
source("../commonSrc/stats/utils/computePercentageExplainedVar.R")
source("../commonSrc/stats/kalmanFilter/computeOneStepAheadObsPredStats.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedMatrices.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedVectors.R")
source("../commonSrc/plot/kalmanFilter/getPlotStateSpectrum.R")
source("../commonSrc/plot/kalmanFilter/getPlotOneStepAheadForecasts.R")
source("../commonSrc/plot/kalmanFilter/getPlotSmoothedStates.R")
source("../commonSrc/plot/kalmanFilter/getPlotPercentageExplainedVar.R")
source("../commonSrc/plot/kalmanFilter/getPlotLogLik.R")

processAll <- function() {
    estNumber <- 26587485
    estMetaDataFilenamePattern <- "results/%08d_estimation.ini"
    estConfigFilenamePattern <- "data/%08d_estimation.ini"
    figFilenamePattern <- "figures/%08d_%s.%s"

    estMetaDataFilename <- sprintf(estMetaDataFilenamePattern, estNumber)
    estMetaData <- read.ini(filepath=estMetaDataFilename)
    estConfigNumber <- as.integer(estMetaData$estimation_config_info$estConfigNumber)

    estConfigFilename <- sprintf(estConfigFilenamePattern, estConfigNumber)
    estConfig <- read.ini(estConfigFilename)
    obsInputMemorySecs <- as.double(estConfig$inputs$obsInputMemorySecs)
    stateInputMemorySecs <- as.double(estConfig$inputs$stateInputMemorySecs)
    stateDim <- eval(parse(text=estConfig$control_variables$stateDim))
    estResFilenamePattern <- estConfig$filenames$estResFilenamePattern

    estResFilename <- sprintf(estResFilenamePattern, estNumber)
    estRes <- get(load(estResFilename))

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

    tSpikeRates <- estRes$tSpikeRates
    time <- estRes$startTime+(1:ncol(estRes$tSpikeRates))/estRes$sRate
    dimObs <- nrow(tSpikeRates)
    nObs <- ncol(tSpikeRates)

    browser()

    show("Plotting logLik")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "logLik", "html")
    fig <- getPlotLogLik(logLik=estRes$kem$iter.record$logLik)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting B")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "B", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "B", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$B, estimated=matrix(coef(kem)$B, nrow=stateDim))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting u")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "u", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "u", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$U, nrow=stateDim))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting Q")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Q", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Q", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting Z")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "Z", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "Z", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$Z, estimated=matrix(coef(kem)$Z, nrow=dimObs))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting a")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "a", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "a", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$A, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting C")
    if(!is.nan(stateInputMemorySecs)) {
        CblockSize <- stateDim*(1+stateInputMemorySamples)
        Coffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Cgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Cgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Cnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Cnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Claser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Claser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CgoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CgoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "CnogoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "CnogoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=stateDim))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
    }

    show("Plotting D")
    if(!is.nan(obsInputMemorySecs)) {
        DblockSize <- dimObs*(1+obsInputMemorySamples)
        Doffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "Dlaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "Dlaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DgoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DgoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, estNumber, "DnogoLaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, "DnogoLaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
    }

    show("Plotting R")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "R", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "R", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$R, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting m0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "m0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "m0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting V0")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "V0", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "V0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$V0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

#     show("Plotting stateSpectrum")
#     pngFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "png")
#     htmlFilename <- sprintf(figFilenamePattern, estNumber, "stateSpectrum", "html")
#     fig <- getPlotStateSpectrum(B=matrix(coef(kem)$B, nrow=stateDim))
#     htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#     # orca(p=fig, file=pngFilename)
#     # print(fig)

    show("Plotting percExpVar")
    percExpVar <- computePercentageExplainedVar(observations=tSpikeRates, predictions=predStats$ytt1)
    pngFilename <- sprintf(figFilenamePattern, estNumber, "percExpVar", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "percExpVar", "html")
    fig <- getPlotPercentageExplainedVar(percExpVar=percExpVar)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting oneStepAheadForecasts")
    if(kem$call$model$R=="diagonal and unequal") {
        R <- diag(coef(kem)$R)
    } else {
        error("Functionality not yet implemented for R!=<diagonal and unequal>")
    }
    predStats <- computeOneStepAheadObsPredStats(xtt1=kfRes$xtt1, Vtt1=kfRes$Vtt1, Z=matrix(coef(kem)$Z, nrow=dimObs), a=as.numeric(coef(kem)$A), D=matrix(coef(kem)$D, nrow=dimObs), R=R, obsInputs=obsInputs)

    pngFilename <- sprintf(figFilenamePattern, estNumber, "oneStepAheadForecasts", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "oneStepAheadForecasts", "html")
    fig <- getPlotOneStepAheadForecasts(time=time, obs=tSpikeRates, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, inputs=obsInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting oneStepAheadForecastsNeuron")
    for(i in 1:nrow(predStats$ytt1)) {
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "html")
        fig <- getPlotOneStepAheadForecasts(time=time, obs=tSpikeRates, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, obsToPlot=c(i), inputs=obsInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting smoothedStates")
    pngFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "png")
    htmlFilename <- sprintf(figFilenamePattern, estNumber, "smoothedStates", "html")
    fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, inputs=stateInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    for(i in 1:nrow(kfRes$xtT)) {
        pngFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smootheState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, estNumber, sprintf("smootheState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, statesToPlot=c(i), inputs=stateInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
    }

    browser()
}

processAll()
