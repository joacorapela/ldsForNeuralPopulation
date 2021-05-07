
require(MARSS)
require(ini)
require(optparse)
require(plotly)
require(R.matlab)
source("../commonSrc/stats/utils/computePercentageExplainedVar.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedMatrices.R")
source("../commonSrc/plot/kalmanFilter/getPlotTrueInitialAndEstimatedVectors.R")
source("../commonSrc/plot/kalmanFilter/getPlotStateSpectrum.R")
source("../commonSrc/plot/kalmanFilter/getStimRecs.R")
source("../commonSrc/plot/kalmanFilter/getPlotSmoothedStates.R")

plotAllRFsAllNeurons <- function(Z, B, C, stateInputMemorySamples, stateInputMemoryToPlotSamples, sRate, mouseName, figFilenamePattern, estNumber, xlab="Delay (sec)", ylab="Value") {
    if(stateInputMemorySamples!=0) {
        stop(sprintf("At the moment only stateInputMemorySamples=0 is supported, but you provided stateInputMemorySamples=%d", stateInputMemorySamples))
    }
    Cvg <- C[,1]
    Cvn <- C[,2]
    Clg <- C[,3]
    Cln <- C[,4]
    ZBCvg <- matrix(NA, nrow=nrow(Z), ncol=stateInputMemoryToPlotSamples+1)
    ZBCvn <- matrix(NA, nrow=nrow(Z), ncol=stateInputMemoryToPlotSamples+1)
    ZBClg <- matrix(NA, nrow=nrow(Z), ncol=stateInputMemoryToPlotSamples+1)
    ZBCln <- matrix(NA, nrow=nrow(Z), ncol=stateInputMemoryToPlotSamples+1)
    Bpow <- diag(rep(1, nrow(B)))
    for(i in 1:(stateInputMemoryToPlotSamples+1)) {
        ZBCvg[,i] <- Z%*%Bpow%*%Cvg
        ZBCvn[,i] <- Z%*%Bpow%*%Cvn
        ZBClg[,i] <- Z%*%Bpow%*%Clg
        ZBCln[,i] <- Z%*%Bpow%*%Cln
        Bpow <- Bpow%*%B
    }

    timeStateInputs <- (0:stateInputMemoryToPlotSamples)/sRate
    nNeurons <- nrow(Z)
    for(n in 1:nNeurons) {
        show(sprintf("Plotting RFs for neuron %d", n))
        fig <- plot_ly(type='scatter', mode='lines+markers')
        fig <- fig%>%add_trace(x=timeStateInputs, y=ZBCvg[n,], name="population_vg", type="scatter", mode="lines+markers")
        fig <- fig%>%add_trace(x=timeStateInputs, y=ZBCvn[n,], name="population_vn", type="scatter", mode="lines+markers")
        fig <- fig%>%add_trace(x=timeStateInputs, y=ZBClg[n,], name="population_lg", type="scatter", mode="lines+markers")
        fig <- fig%>%add_trace(x=timeStateInputs, y=ZBCln[n,], name="population_ln", type="scatter", mode="lines+markers")
        fig <- fig %>% layout(xaxis=list(title=xlab), yaxis=list(title=ylab))
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("allRFsNeuron%d", n), "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("allRFsNeuron%d", n), "html")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
    }
}

processAll <- function() {
    # DSSSM
    # x_n = B x_{n-1} + u + C c_n + G w_n, w_n ~ N(0, Q)
    # y_n = Z x_n     + a + D d_n + H v_n, v_n ~ N(0, R)
    
    # PLDS
    # x_n = A x_{n-1} + B u_n + w_n, w_n ~ N(0, Q)
    # z_n = C x_n     + D d_n + d
    # y_{ni} ~ P(\lambda=exp(z_{ni})
   
    mouseName <- "MPV18_2"
    estNumber <- 0
    startTimeSecs <- 180
    stateInputMemorySecs <- 0.0
    stateInputMemoryToPlotSamples <- 10
    timeSeriesFilenamePattern <- "../../data/%s/binLDStimeSeries.ini"
    timeSeriesFilename <- sprintf(timeSeriesFilenamePattern, mouseName)
    timeSeriesConfig <- read.ini(timeSeriesFilename)
    dataFilename <- timeSeriesConfig$filenames$saveFilename
    resultsFilename <- "../../results/MPV18_2/PLDSresults.mat"
    figFilenamePattern <- "../../figures/%s/%08d_%s.%s"

    timeSeries <- readMat(dataFilename)$timeSeries

    sRate <- timeSeries["sRate", 1, 1][[1]][1, 1]

    loadRes <- readMat(resultsFilename)

    dsSSM <- list()
    dsSSM$B <- loadRes$params["model",1,1][[1]]["A",1,1][[1]]
    dsSSM$C <- loadRes$params["model",1,1][[1]]["B",1,1][[1]]
    dsSSM$Q <- loadRes$params["model",1,1][[1]]["Q",1,1][[1]]
    dsSSM$m0 <- loadRes$params["model",1,1][[1]]["x0",1,1][[1]]
    dsSSM$V0 <- loadRes$params["model",1,1][[1]]["Q0",1,1][[1]]
    dsSSM$Z <- loadRes$params["model",1,1][[1]]["C",1,1][[1]]
    dsSSM$a <- loadRes$params["model",1,1][[1]]["d",1,1][[1]]
    trainSpikeCounts <- loadRes$y

    spikeCounts <- loadRes$y
    stateInputs <- loadRes$u
    stateDim <- nrow(dsSSM$B)
    obsDim <- nrow(spikeCounts)
    nObs <- ncol(spikeCounts)
    stateInputMemorySamples <- stateInputMemorySecs*sRate
    u <- matrix(rep(0, stateDim), ncol=1)
    # D <- matrix(rep(0, obsDim), ncol=1)
    # d <- matrix(rep(0, nObs), nrow=1)

    xnN <- loadRes$seq["posterior",,][[1]]["xsm",1,1][[1]]
    dim(xnN) <- c(nrow(xnN), 1, ncol(xnN))
    VnN <- array(loadRes$seq["posterior",,][[1]]["Vsm",1,1][[1]], dim=c(stateDim, stateDim, nObs))
    ksRes <- list(xnN=xnN, VnN=VnN)

    time <- startTimeSecs+(1:ncol(spikeCounts))/sRate
    minTime <- min(time)
    maxTime <- max(time)

    goStimOnSecs <- timeSeries["goStimOnSecs", 1, 1][[1]][1,]
    goStimOffSecs <- timeSeries["goStimOffSecs", 1, 1][[1]][1,]
    toKeepIndices <- which(minTime<=goStimOnSecs & goStimOffSecs<=maxTime)
    goStimOnSecs <- goStimOnSecs[toKeepIndices]
    goStimOffSecs <- goStimOffSecs[toKeepIndices]

    nogoStimOnSecs <- timeSeries["nogoStimOnSecs", 1, 1][[1]][1,]
    nogoStimOffSecs <- timeSeries["nogoStimOffSecs", 1, 1][[1]][1,]
    toKeepIndices <- which(minTime<=nogoStimOnSecs & nogoStimOffSecs<=maxTime)
    nogoStimOnSecs <- nogoStimOnSecs[toKeepIndices]
    nogoStimOffSecs <- nogoStimOffSecs[toKeepIndices]

    laserStimOnSecs <- timeSeries["laserStimOnSecs", 1, 1][[1]][1,]
    laserStimOffSecs <- timeSeries["laserStimOffSecs", 1, 1][[1]][1,]
    toKeepIndices <- which(minTime<=laserStimOnSecs & laserStimOffSecs<=maxTime)
    laserStimOnSecs <- laserStimOnSecs[toKeepIndices]
    laserStimOffSecs <- laserStimOffSecs[toKeepIndices]

if(FALSE) {
    show("Plotting logLik")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "logLik", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "logLik", "html")
    fig <- getPlotLogLik(logLik=dsSSM$logLik)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)
}

if(FALSE) {
    show("Plotting B")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "B", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "B", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=dsSSM$B)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting stateSpectrum")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "stateSpectrum", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "stateSpectrum", "html")
    fig <- getPlotStateSpectrum(B=dsSSM$B)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

#     show("Plotting U")
#     pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "U", "png")
#     htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "U", "html")
#     fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$u)
#     htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#     orca(p=fig, file=pngFilename)
    # print(fig)

    if(!is.nan(stateInputMemorySecs)) {
        show("Plotting C")
        x <- (0:stateInputMemorySamples)/sRate
        CblockSize <- 1+stateInputMemorySamples
        Coffset <- 0
        #
        CVisualStimGo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CVisualStimNogo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CVisualStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimGo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimGo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Coffset <- Coffset + CblockSize
        #
        CLaserStimNogo <- matrix(data=dsSSM$C[,Coffset+(1:CblockSize)], ncol=CblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "CLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=CLaserStimNogo, estimatedLegendLabelPattern="state %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

    show("Plotting Q")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "Q", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "Q", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=diag(dsSSM$Q), xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting m0")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "m0", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "m0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$m0, xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting V0")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "V0", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "V0", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=diag(dsSSM$V0), xlab="State Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting Z")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "Z", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "Z", "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(estimated=dsSSM$Z)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    show("Plotting a")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "a", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "a", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=dsSSM$a, xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)
}

if(FALSE) {
    if(!is.nan(obsInputMemorySecs)) {
        show("Plotting D")
        x <- (0:obsInputMemorySamples)/sRate
        DblockSize <- 1+obsInputMemorySamples
        Doffset <- 0
        #
        DVisualStimGo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DVisualStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DVisualStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DVisualStimGo, estimatedLegendLabelPattern="neuron %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        DVisualStimNogo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DVisualStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DVisualStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DVisualStimNogo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset + DblockSize
        #
        DLaserStimGo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DLaserStimGo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DLaserStimGo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DLaserStimGo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        Doffset <- Doffset+DblockSize
        #
        DLaserStimNogo <- matrix(dsSSM$D[,Doffset+(1:DblockSize)], ncol=DblockSize)
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DLaserStimNogo", "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "DLaserStimNogo", "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(x=x, estimated=DLaserStimNogo, estimatedLegendLabelPattern="observation %d", xlab="Delay (sec)")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }
}
if(FALSE) {
    show("Plotting R")
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "R", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "R", "html")
    fig <- getPlotTrueInitialAndEstimatedVectors(estimated=diag(dsSSM$R), xlab="Observation Index")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

    predStats <- computeOneStepAheadObsPredStats(xtt1=kfRes$xnn1[,1,], Vtt1=kfRes$Vnn1, Z=dsSSM$Z, a=as.vector(dsSSM$a), D=dsSSM$D, R=dsSSM$R, obsInputs=estRes$obsInputs[,1,])

    show("Plotting percExpVar")
    percExpVar <- computePercentageExplainedVar(observations=trainSqrtSpikeCounts, predictions=predStats$ytt1)
    pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "percExpVar", "png")
    htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "percExpVar", "html")
    fig <- getPlotPercentageExplainedVar(percExpVar=percExpVar)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    orca(p=fig, file=pngFilename)
    # print(fig)

#     show("Plotting oneStepAheadForecasts")
#     pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "oneStepAheadForecasts", "png")
#     htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "oneStepAheadForecasts", "html")
#     fig <- getPlotOneStepAheadForecasts(time=time, obs=trainSqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff)
#     htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#     orca(p=fig, file=pngFilename)
    # print(fig)

#     for(i in 1:nrow(predStats$ytt1)) {
#         show(sprintf("Plotting oneStepAheadForecast for neuron %d", i))
#         pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "png")
#         htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("oneStepAheadForecastsNeuron%d", i), "html")
#         fig <- getPlotOneStepAheadForecasts(time=time, obs=trainSqrtSpikeCounts, ytt1=predStats$ytt1, Wtt1=predStats$Wtt1, goStimOn=goStimOn, goStimOff=goStimOff, nogoStimOn=nogoStimOn, nogoStimOff=nogoStimOff, laserStimOn=laserStimOn, laserStimOff=laserStimOff, obsToPlot=c(i))
#         htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#         orca(p=fig, file=pngFilename)
        # print(fig)
#     }
# }
}
#     show("Plotting smoothedStates")
#     pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "smoothedStates", "png")
#     htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, "smoothedStates", "html")
#     fig <- getPlotSmoothedStates(time=time, xtT=ksRes$xnN[,1,], VtT=ksRes$VnN, goStimOn=goStimOnSecs, goStimOff=goStimOffSecs, nogoStimOn=nogoStimOnSecs, nogoStimOff=nogoStimOffSecs, laserStimOn=laserStimOnSecs, laserStimOff=laserStimOffSecs)
#     htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
#     orca(p=fig, file=pngFilename)
#     # print(fig)

    for(i in 1:nrow(ksRes$xnN)) {
        show(sprintf("Plotting smoothedState %d", i))
        pngFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("smoothedState%d", i), "png")
        htmlFilename <- sprintf(figFilenamePattern, mouseName, estNumber, sprintf("smoothedState%d", i), "html")
        fig <- getPlotSmoothedStates(time=time, xtT=ksRes$xnN[,1,], VtT=ksRes$VnN, goStimOn=goStimOnSecs, goStimOff=goStimOffSecs, nogoStimOn=nogoStimOnSecs, nogoStimOff=nogoStimOffSecs, laserStimOn=laserStimOnSecs, laserStimOff=laserStimOffSecs, statesToPlot=c(i))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
    }

#     plotAllRFsAllNeurons(Z=dsSSM$Z, B=dsSSM$B, C=dsSSM$C, stateInputMemorySamples=stateInputMemorySamples, stateInputMemoryToPlotSamples=stateInputMemoryToPlotSamples, sRate=sRate, figFilenamePattern=figFilenamePattern, mouseName=mouseName, estNumber=estNumber)

    browser()
}

processAll()
