require(R.matlab)
require(plotly)                                                    

getRasterPlot <- function(spikesTimes, stimOnTimeStamps, stimOffTimeStamps, goTrialIndices, nogoTrialIndices, laserOnsets, laserDuration, sRate, startTime, endTime, markerColor="rgb(128,128,128)", goStimColor="rgb(0,255,0)", nogoStimColor="rgb(255,0,0)", laserColor="rgb(0,0,255)", markerSize=3, stimOpacity=0.3, laserOpacity=0.3) {

    getStimRecs <- function(stimOnTimeStamps, stimOffTimeStamps, trialIndices, sRate, color, opacity, ymin, ymax) {
        recs <- vector(mode="list", length=length(trialIndices))
        for(i in 1:length(trialIndices)) {
            x0 <- stimOnTimeStamps[trialIndices[i]]/sRate
            x1 <- stimOffTimeStamps[trialIndices[i]]/sRate
            recs[[i]] <- list(type="rect",
                                    fillcolor=color, 
                                    line=list(color=color), 
                                    opacity=opacity,
                                    x0=x0, x1=x1, xref="x",
                                    y0=ymin, y1=ymax, yref="y")
        }
        return(recs)
    }
    getLaserRecs <- function(laserOnsets, laserDuration, color, opacity, ymin, ymax) {
        recs <- vector(mode="list", length=length(laserOnsets))
        for(i in 1:length(laserOnsets)) {
            x0 <- laserOnsets[i]
            x1 <- laserOnsets[i]+laserDuration
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
    goRecs <- getStimRecs(stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, trialIndices=goValidTrialIndices, sRate=sRate, color=goStimColor, opacity=stimOpacity, ymin=0, ymax=nUnits)
    nogoValidTrialIndices <- getValidTrialIndices(stimTrialIndices=nogoTrialIndices, stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, startTime=startTime, endTime=endTime, sRate=sRate)
    nogoRecs <- getStimRecs(stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, trialIndices=nogoValidTrialIndices, sRate=sRate, color=nogoStimColor, opacity=stimOpacity, ymin=0, ymax=nUnits)
   
    validLaserOnsets <- laserOnsets[!is.nan(laserOnsets) & startTime<=laserOnsets & laserOnsets+0.1<endTime]
    laserRecs <- getLaserRecs(laserOnsets=validLaserOnsets, laserDuration=laserDuration, color=laserColor, opacity=laserOpacity, ymin=0, ymax=nUnits)
    stimShapes <- c(goRecs, nogoRecs, laserRecs)

    fig <- fig%>%layout(shapes=stimShapes, xaxis=list(title="Time (sec)"), yaxis=list(title="Neuron Index"))

    return(fig)
}

processAll <- function() {
    sRate <- 30000
    v1ShaftIndex <- 1
    lmShaftIndex <- 1
    laserDuration <- .1
    matlabDataFilename <- "../../data/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_subset_V6.mat"
    figFilenamePattern <- "figures/%sShaft%dRasterPlot-start%.02fsec-end%.02fsec.%s"

    loadRes <- readMat(matlabDataFilename)


    aux  <- loadRes[[sprintf("V1Shaft%d", v1ShaftIndex)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    v1UnitsSpikesTimeStamps <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        v1UnitsSpikesTimeStamps[[n]] <- aux[[n]][[1]][,1]
    }

    aux  <- loadRes[[sprintf("LMShaft%d", v1ShaftIndex)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    lmUnitsSpikesTimeStamps <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        lmUnitsSpikesTimeStamps[[n]] <- aux[[n]][[1]][,1]
    }

    stimOnTimeStamps <- loadRes[["PStepTimeStampOn"]][1,]
    stimOffTimeStamps <- loadRes[["PStepTimeStampOff"]][1,]
    goTrialIndices <- loadRes[["gotrialind"]][,1]
    nogoTrialIndices <- loadRes[["nogotrialind"]][,1]
    laserDelays <- loadRes[["LaserDelay"]][1,]/1000 #laser delays are in ms and I want them in sec

    laserOnsets <- stimOnTimeStamps/sRate+laserDelays
    # startTimeToPlot <- 3000
    startTimeToPlot <- 0
    plotDuration <- 240
    endTimeToPlot <- startTimeToPlot+plotDuration

    v1Fig <- getRasterPlot(spikesTimes=v1UnitsSpikesTimeStamps, stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, goTrialIndices=goTrialIndices, nogoTrialIndices=nogoTrialIndices, laserOnsets=laserOnsets, laserDuration=laserDuration, sRate=sRate, startTime=startTimeToPlot, endTime=endTimeToPlot)
    # pngFigFilename <- sprintf(figFilenamePattern, "v1", v1ShaftIndex, startTimeToPlot, endTimeToPlot, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, "v1", v1ShaftIndex, startTimeToPlot, endTimeToPlot, "html")
    # orca(p=v1Fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(v1Fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(v1Fig)

    lmFig <- getRasterPlot(spikesTimes=lmUnitsSpikesTimeStamps, stimOnTimeStamps=stimOnTimeStamps, stimOffTimeStamps=stimOffTimeStamps, goTrialIndices=goTrialIndices, nogoTrialIndices=nogoTrialIndices, laserOnsets=laserOnsets, laserDuration=laserDuration, sRate=sRate, startTime=startTimeToPlot, endTime=endTimeToPlot)
    # pngFigFilename <- sprintf(figFilenamePattern, "lm", lmShaftIndex, startTimeToPlot, endTimeToPlot, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, "lm", lmShaftIndex, startTimeToPlot, endTimeToPlot, "html")
    # orca(p=lmFig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(lmFig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(lmFig)

    browser()
}

processAll()
rm(processAll)
