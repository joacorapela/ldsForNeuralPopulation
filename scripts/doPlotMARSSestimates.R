
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
    dimState <- 3
    stateInputMemorySecs <- 0
    obsInputMemorySecs <- 0
    # resultsFilename <- "results/dimState03_stateInputMemory0.00_obsInputMemory0.00_kas.RData"
    # figFilenamePattern <- "figures/dimState%02d_stateInputMemory%.02f_obsInputMemory%.02f_kas_%s.%s"
    resultsFilename <- "results/dimState03_stateInputMemory0.00_obsInputMemory0.00_ss.RData"
    figFilenamePattern <- "figures/dimState%02d_stateInputMemory%.02f_obsInputMemory%.02f_ss_%s.%s"

    results <- get(load(resultsFilename))
    kem <- results$kem
    obsInputs <- results$obsInputs
    stateInputs <- results$stateInputs
    initialConds <- list(B=matrix(kem$start$B, ncol=dimState), Z=matrix(kem$start$Z, ncol=dimState), RDiag=as.vector(kem$start$R))
    sRate <- results$sRate
    stateInputMemorySamples <- stateInputMemorySecs*sRate
    obsInputMemorySamples <- obsInputMemorySecs*sRate

    kfRes <- MARSSkf(kem)
#     kem0 <- kem
#     kem0$par <- kem0$start
#     kfRes0 <- MARSSkf(kem0)

    tSpikeRates <- results$tSpikeRates
    time <- results$startTime+(1:ncol(results$tSpikeRates))/results$sRate
    dimObs <- nrow(tSpikeRates)
    nObs <- ncol(tSpikeRates)

    browser()

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "logLik", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "logLik", "html")
    fig <- getPlotLogLik(logLik=results$kem$iter.record$logLik)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "B", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "B", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$B, estimated=matrix(coef(kem)$B, nrow=dimState))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "u", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "u", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$U, nrow=dimState))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Q", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Q", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Z", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Z", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$Z, estimated=matrix(coef(kem)$Z, nrow=dimObs))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "a", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "a", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$A, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(stateInputMemorySecs)) {
        CblockSize <- dimState*(1+stateInputMemorySamples)
        Coffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Cgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Cgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=dimState))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Cnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Cnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=dimState))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Claser", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Claser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$C[Coffset+(1:CblockSize)], nrow=dimState))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
    }

    if(!is.nan(obsInputMemorySecs)) {
        DblockSize <- dimObs*(1+obsInputMemorySamples)
        Doffset <- 0
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dgo", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dgo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dnogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dnogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dlaser", "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "Dlaser", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=matrix(coef(kem)$D[Doffset+(1:DblockSize)], nrow=dimObs))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
    }

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "R", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "R", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$R, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "m0", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "m0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "V0", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "V0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$V0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

#     pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "stateSpectrum", "png")
#     htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "stateSpectrum", "html")
#     fig <- getPlotStateSpectrum(B=matrix(coef(kem)$B, nrow=dimState))
#     htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#     # orca(p=fig, file=pngFilename)
#     # print(fig)

    if(kem$call$model$R=="diagonal and unequal") {
        R <- diag(coef(kem)$R)
    } else {
        error("Functionality not yet implemented for R!=<diagonal and unequal>")
    }
    predStats <- computeOneStepAheadObsPredStats(xtt1=kfRes$xtt1, Vtt1=kfRes$Vtt1, Z=matrix(coef(kem)$Z, nrow=dimObs), a=as.numeric(coef(kem)$A), D=matrix(coef(kem)$D, nrow=dimObs), R=R, obsInputs=obsInputs)

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "oneStepAheadForecasts", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "oneStepAheadForecasts", "html")
    fig <- getPlotOneStepAheadForecasts(time=time, obs=tSpikeRates, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, inputs=obsInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    percExpVar <- computePercentageExplainedVar(observations=tSpikeRates, predictions=predStats$ytt1)
    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "percExpVar", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "percExpVar", "html")
    fig <- getPlotPercentageExplainedVar(percExpVar=percExpVar)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)


    for(i in 1:nrow(predStats$ytt1)) {
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, sprintf("oneStepAheadForecastsNeuron%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, sprintf("oneStepAheadForecastsNeuron%d", i), "html")
        fig <- getPlotOneStepAheadForecasts(time=time, obs=tSpikeRates, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, obsToPlot=c(i), inputs=obsInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
    }

    pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "smoothedStates", "png")
    htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, "smoothedStates", "html")
    fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, inputs=stateInputs)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # orca(p=fig, file=pngFilename)
    # print(fig)

    for(i in 1:nrow(kfRes$xtT)) {
        pngFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, sprintf("smootheState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, dimState, stateInputMemorySecs, obsInputMemorySecs, sprintf("smootheState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=kfRes$xtT, VtT=kfRes$VtT, statesToPlot=c(i), inputs=stateInputs)
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        # orca(p=fig, file=pngFilename)
        # print(fig)
    }

    browser()
}

processAll()
