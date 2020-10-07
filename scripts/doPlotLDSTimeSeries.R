require(R.matlab)
require(plotly)

getRasterPlot <- function(spikesTimes, stimOnTimeStamps, stimOffTimeStamps, goTrialIndices, nogoTrialIndices, laserOnsets, laserDuration, sRate, startTime, endTime, markerColor="rgb(128,128,128)", goStimColor="rgb(0,255,0)", nogoStimColor="rgb(255,0,0)", laserColor="rgb(0,0,255)", markerSize=3, stimOpacity=0.3, laserOpacity=0.3) {

    getRecs <- function(onsets, durations, color, opacity, ymin, ymax) {
        recs <- vector(mode="list", length=length(laserOnsets))
        for(i in 1:length(onsets)) {
            x0 <- onsets[i]
            x1 <- onsets[i]+durations[i]
            recs[[i]] <- list(type="rect",
                              fillcolor=color,
                              line=list(color=color),
                              opacity=opacity,
                              x0=x0, x1=x1, xref="x",
                              y0=ymin, y1=ymax, yref="y")
        }
        return(recs)
    }

    getValidTrialIndices <- function(stimTrialIndices, stimOnTimeStamps, stimOffTimeStamps, startTime, endTime, sRate) {
        validStimTrialIndices <- c()
        for(trialIndex in stimTrialIndices) {
            if(!(stimOffTimeStamps[trialIndex]/sRate<startTime || stimOnTimeStamps[trialIndex]/sRate>endTime)) {
                validStimTrialIndices <- c(validStimTrialIndices, trialIndex)
            }
        }
        return(validStimTrialIndices)
    }

    nUnits <- length(spikesTimes)

    # raster plot
    fig <- plot_ly(type="scatter", mode="markers")
    for(n in 1:nUnits) {
        unitSpikeTimesSec <- spikesTimes[[n]]/sRate
        spikeTimesToPlot <- unitSpikeTimesSec[startTime<=unitSpikeTimesSec & unitSpikeTimesSec<=endTime]
        traceName <- sprintf("n%02d", n)
        fig <- fig%>%add_trace(x=spikeTimesToPlot, y=rep(n, times=length(spikeTimesToPlot)), marker=list(color=markerColor, size=markerSize), name=traceName, showlegend=TRUE)
    }

    # stim shapes
    goValidTrialIndices <- getValidTrialIndices(stimTrialIndices=goTrialIndices, stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, startTime=startTime, endTime=endTime, sRate=sRate)
    goRecs <- getRecs(onsets=stimOnTimeStamps[goValidTrialIndices]/sRate, durations=(stimOffTimeStamps[goValidTrialIndices]-stimOnTimeStamps[goValidTrialIndices])/sRate, color=goStimColor, opacity=stimOpacity, ymin=0, ymax=nUnits)
    nogoValidTrialIndices <- getValidTrialIndices(stimTrialIndices=nogoTrialIndices, stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, startTime=startTime, endTime=endTime, sRate=sRate)
    nogoRecs <- getRecs(onsets=stimOnTimeStamps[nogoValidTrialIndices]/sRate, durations=(stimOffTimeStamps[nogoValidTrialIndices]-stimOnTimeStamps[nogoValidTrialIndices])/sRate, color=nogoStimColor, opacity=stimOpacity, ymin=0, ymax=nUnits)

    validLaserOnsets <- laserOnsets[!is.nan(laserOnsets) & startTime<=laserOnsets & laserOnsets+0.1<endTime]
    laserRecs <- getRecs(onsets=validLaserOnsets, durations=rep(laserDuration, times=length(validLaserOnsets)), color=laserColor, opacity=laserOpacity, ymin=0, ymax=nUnits)
    stimShapes <- c(goRecs, nogoRecs, laserRecs)

    fig <- fig%>%layout(shapes=stimShapes, xaxis=list(title="Time (sec)"), yaxis=list(title="Neuron Index"))

    return(fig)
}

processAll <- function() {
    visualGoColor <-     "rgb(0,   255, 0)"
    visualNogoColor <-   "rgb(255, 0,   0)"
    laserGoColor <-  "rgb(0,   0,   255)"
    laserNogoColor <-  "rgb(0,   255,   255)"
    spikesColor <- "rgb(0,   0,   0)"
    v1Shaft1NeuronToPlot <- 22
    v1Shaft2NeuronToPlot <- 16
    lmShaft1NeuronToPlot <- 5
    lmShaft2NeuronToPlot <- 4
    sRate <- 30000
    startTimeToPlot <- 180
    # startTimeToPlot <- 50*60
    plotDuration <- 180
    endTimeToPlot <- startTimeToPlot+plotDuration
    timeSeriesFilename <- "results/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_timeSeries.RData"
    figFilenamePattern <- "figures/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_%sShaft%dNeuron%d_start%.02fsec_end%.02fsec.%s"

    loadRes <- get(load(timeSeriesFilename))
    breaks <- loadRes$breaks
    times <- (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2/sRate
    samplesToPlot <- which(startTimeToPlot<=times & times<endTimeToPlot)

    goStim <- loadRes$goStim
    figGoStim <- plot_ly(x=times[samplesToPlot], y=goStim[samplesToPlot], name="visual go", line=list(color=visualGoColor), marker=list(color=visualGoColor), type="scatter", mode="lines+markers")
    figGoStim <- figGoStim%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Visual Go Stimulus"))

    nogoStim <- loadRes$nogoStim
    figNogoStim <- plot_ly(x=times[samplesToPlot], y=nogoStim[samplesToPlot], name="visual nogo", line=list(color=visualNogoColor), marker=list(color=visualNogoColor), type="scatter", mode="lines+markers")
    figNogoStim <- figNogoStim%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Visual Nogo Stimulus"))

    laserStim <- loadRes$laserStim
    laserGoStim <- laserStim*goStim
    laserNogoStim <- laserStim*nogoStim
    figLaserGo <- plot_ly(x=times[samplesToPlot], y=laserGoStim[samplesToPlot], name="laser go", line=list(color=laserGoColor), marker=list(color=laserGoColor), type="scatter", mode="lines+markers")
    figLaserGo <- figLaserGo%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Laser Go Stimulus"))
    figLaserNogo <- plot_ly(x=times[samplesToPlot], y=laserNogoStim[samplesToPlot], name="laser nogo", line=list(color=laserNogoColor), marker=list(color=laserNogoColor), type="scatter", mode="lines+markers")
    figLaserNogo <- figLaserNogo%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Laser Nogo Stimulus"))

    v1Shaft1SpikeCounts <- loadRes$v1Shaft1SpikeCounts
    v1Shaft1FigSpikeCounts <- plot_ly(x=times[samplesToPlot], y=v1Shaft1SpikeCounts[v1Shaft1NeuronToPlot, samplesToPlot], marker=list(color=spikesColor), name="spikes", type="scatter", mode="markers")
    v1Shaft1FigSpikeCounts <- v1Shaft1FigSpikeCounts%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Spike Count"))
    v1Shaft1AllPlots <- list(figLaserGo, figLaserNogo, figGoStim, figNogoStim, v1Shaft1FigSpikeCounts)
    v1Shaft1Fig <- subplot(v1Shaft1AllPlots, nrows=length(v1Shaft1AllPlots), shareX=TRUE)
    #
    v1Shaft1PNGFigFilename <- sprintf(figFilenamePattern, "v1TimeSeries", 1, v1Shaft1NeuronToPlot, startTimeToPlot, endTimeToPlot, "png")
    v1Shaft1HTMLFigFilename <- sprintf(figFilenamePattern, "v1TimeSeries", 1, v1Shaft1NeuronToPlot, startTimeToPlot, endTimeToPlot, "html")
    orca(p=v1Shaft1Fig, file=v1Shaft1PNGFigFilename)
    htmlwidgets::saveWidget(as_widget(v1Shaft1Fig), file.path(normalizePath(dirname(v1Shaft1HTMLFigFilename)), basename(v1Shaft1HTMLFigFilename)))
    # print(v1Shaft1Fig)

    v1Shaft2SpikeCounts <- loadRes$v1Shaft2SpikeCounts
    v1Shaft2FigSpikeCounts <- plot_ly(x=times[samplesToPlot], y=v1Shaft2SpikeCounts[v1Shaft2NeuronToPlot, samplesToPlot], marker=list(color=spikesColor), name="spikes", type="scatter", mode="markers")
    v1Shaft2FigSpikeCounts <- v1Shaft2FigSpikeCounts%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Spike Count"))
    v1Shaft2AllPlots <- list(figLaserGo, figLaserNogo, figGoStim, figNogoStim, v1Shaft2FigSpikeCounts)
    v1Shaft2Fig <- subplot(v1Shaft2AllPlots, nrows=length(v1Shaft2AllPlots), shareX=TRUE)
    #
    v1Shaft2PNGFigFilename <- sprintf(figFilenamePattern, "v1TimeSeries", 2, v1Shaft2NeuronToPlot, startTimeToPlot, endTimeToPlot, "png")
    v1Shaft2HTMLFigFilename <- sprintf(figFilenamePattern, "v1TimeSeries", 2, v1Shaft2NeuronToPlot, startTimeToPlot, endTimeToPlot, "html")
    orca(p=v1Shaft2Fig, file=v1Shaft2PNGFigFilename)
    htmlwidgets::saveWidget(as_widget(v1Shaft2Fig), file.path(normalizePath(dirname(v1Shaft2HTMLFigFilename)), basename(v1Shaft2HTMLFigFilename)))
    # print(v1Shaft2Fig)

    lmShaft1SpikeCounts <- loadRes$lmShaft1SpikeCounts
    lmShaft1FigSpikeCounts <- plot_ly(x=times[samplesToPlot], y=lmShaft1SpikeCounts[lmShaft1NeuronToPlot, samplesToPlot], marker=list(color=spikesColor), name="spikes", type="scatter", mode="markers")
    lmShaft1FigSpikeCounts <- lmShaft1FigSpikeCounts%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Spike Count"))
    lmShaft1AllPlots <- list(figLaserGo, figLaserNogo, figGoStim, figNogoStim, lmShaft1FigSpikeCounts)
    lmShaft1Fig <- subplot(lmShaft1AllPlots, nrows=length(lmShaft1AllPlots), shareX=TRUE)
    #
    lmShaft1PNGFigFilename <- sprintf(figFilenamePattern, "lmTimeSeries", 1, lmShaft1NeuronToPlot, startTimeToPlot, endTimeToPlot, "png")
    lmShaft1HTMLFigFilename <- sprintf(figFilenamePattern, "lmTimeSeries", 1, lmShaft1NeuronToPlot, startTimeToPlot, endTimeToPlot, "html")
    orca(p=lmShaft1Fig, file=lmShaft1PNGFigFilename)
    htmlwidgets::saveWidget(as_widget(lmShaft1Fig), file.path(normalizePath(dirname(lmShaft1HTMLFigFilename)), basename(lmShaft1HTMLFigFilename)))
    # print(lmShaft1Fig)

    lmShaft2SpikeCounts <- loadRes$lmShaft2SpikeCounts
    lmShaft2FigSpikeCounts <- plot_ly(x=times[samplesToPlot], y=lmShaft2SpikeCounts[lmShaft2NeuronToPlot, samplesToPlot], marker=list(color=spikesColor), name="spikes", type="scatter", mode="markers")
    lmShaft2FigSpikeCounts <- lmShaft2FigSpikeCounts%>%layout(xaxis=list(title="Time (sec)"), yaxis=list(title="Spike Count"))
    lmShaft2AllPlots <- list(figLaserGo, figLaserNogo, figGoStim, figNogoStim, lmShaft2FigSpikeCounts)
    lmShaft2Fig <- subplot(lmShaft2AllPlots, nrows=length(lmShaft2AllPlots), shareX=TRUE)
    #
    lmShaft2PNGFigFilename <- sprintf(figFilenamePattern, "lmTimeSeries", 2, lmShaft2NeuronToPlot, startTimeToPlot, endTimeToPlot, "png")
    lmShaft2HTMLFigFilename <- sprintf(figFilenamePattern, "lmTimeSeries", 2, lmShaft2NeuronToPlot, startTimeToPlot, endTimeToPlot, "html")
    orca(p=lmShaft2Fig, file=lmShaft2PNGFigFilename)
    htmlwidgets::saveWidget(as_widget(lmShaft2Fig), file.path(normalizePath(dirname(lmShaft2HTMLFigFilename)), basename(lmShaft2HTMLFigFilename)))
    # print(lmShaft2Fig)

    browser()
}

processAll()
rm(processAll)
