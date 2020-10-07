require(R.matlab)
require(plotly)                                                    

getSpikeCounts <- function(spikeSamples, breaks) {
    nUnits <- length(spikeSamples)
    spikeCounts <- matrix(NA, nrow=nUnits, ncol=length(breaks)-1)
    for(n in 1:nUnits) {
        spikeCounts[n,] <- hist(spikeSamples[[n]], breaks=breaks, plot=FALSE)$counts
    }
    return(spikeCounts)
}

getBinnedStimulus <- function(stimOnSamples, stimOffSamples, breaks) {
    binnedStimulus <- rep(0, times=length(breaks)-1)
    for(i in 1:length(stimOnSamples)) {
        stimulatedBins <- which(stimOnSamples[i]<=breaks & breaks<stimOffSamples[i])
        binnedStimulus[stimulatedBins] <- 1.0
    }
    return(binnedStimulus)
}

processAll <- function() {
    sRate <- 30000
    binSizeSecs <- 0.1 # in sec
    laserDuration <- .1
    matlabDataFilename <- "../../data/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_subset_V6.mat"
    saveFilename <- "results/task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_timeSeries.RData"

    loadRes <- readMat(matlabDataFilename)

    maxSpikeSample <- 0

    aux  <- loadRes[[sprintf("V1Shaft%d", 1)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    v1Shaft1SpikeSamples <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        v1Shaft1SpikeSamples[[n]] <- aux[[n]][[1]][,1]
        maxSpikeSample <- max(maxSpikeSample, max(v1Shaft1SpikeSamples[[n]]))
    }

    aux  <- loadRes[[sprintf("V1Shaft%d", 2)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    v1Shaft2SpikeSamples <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        v1Shaft2SpikeSamples[[n]] <- aux[[n]][[1]][,1]
        maxSpikeSample <- max(maxSpikeSample, max(v1Shaft2SpikeSamples[[n]]))
    }

    aux  <- loadRes[[sprintf("LMShaft%d", 1)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    lmShaft1SpikeSamples <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        lmShaft1SpikeSamples[[n]] <- aux[[n]][[1]][,1]
        maxSpikeSample <- max(maxSpikeSample, max(lmShaft1SpikeSamples[[n]]))
    }

    aux  <- loadRes[[sprintf("LMShaft%d", 2)]]["SingleUnitSpikeTimes",1,1][[1]]
    nUnits <- length(aux)
    lmShaft2SpikeSamples <- vector(mode="list", length=nUnits)
    for(n in 1:nUnits) {
        lmShaft2SpikeSamples[[n]] <- aux[[n]][[1]][,1]
        maxSpikeSample <- max(maxSpikeSample, max(lmShaft2SpikeSamples[[n]]))
    }

    binSizeSamples <- round(binSizeSecs*sRate)
    breaks <- seq(from=0, to=maxSpikeSample+binSizeSamples, by=binSizeSamples)
    v1Shaft1SpikeCounts <- getSpikeCounts(spikeSamples=v1Shaft1SpikeSamples, breaks=breaks)
    v1Shaft2SpikeCounts <- getSpikeCounts(spikeSamples=v1Shaft2SpikeSamples, breaks=breaks)
    lmShaft1SpikeCounts <- getSpikeCounts(spikeSamples=lmShaft1SpikeSamples, breaks=breaks)
    lmShaft2SpikeCounts <- getSpikeCounts(spikeSamples=lmShaft2SpikeSamples, breaks=breaks)

    stimOnSamples <- loadRes[["PStepTimeStampOn"]][1,]
    stimOffSamples <- loadRes[["PStepTimeStampOff"]][1,]
    goTrialIndices <- loadRes[["gotrialind"]][,1]
    nogoTrialIndices <- loadRes[["nogotrialind"]][,1]
    laserDelays <- loadRes[["LaserDelay"]][1,]/1000 #laser delays are in ms and I want them in sec

    goStimBinned <- getBinnedStimulus(stimOnSamples=stimOnSamples[goTrialIndices], stimOffSamples=stimOffSamples[goTrialIndices], breaks=breaks)
    nogoStimBinned <- getBinnedStimulus(stimOnSamples=stimOnSamples[nogoTrialIndices], stimOffSamples=stimOffSamples[nogoTrialIndices], breaks=breaks)
    laserOnsetSamples <- stimOnSamples+round(laserDelays*sRate)
    laserOnsetSamples <- laserOnsetSamples[which(!is.nan(laserOnsetSamples))]
    laserStimBinned <- getBinnedStimulus(stimOnSamples=laserOnsetSamples, stimOffSamples=laserOnsetSamples+round(laserDuration*sRate), breaks=breaks)

    timeSeries <- list(sRate=1.0/binSizeSecs, breaks=breaks, v1Shaft1SpikeCounts=v1Shaft1SpikeCounts, v1Shaft2SpikeCounts=v1Shaft2SpikeCounts, lmShaft1SpikeCounts=lmShaft1SpikeCounts, lmShaft2SpikeCounts=lmShaft2SpikeCounts, goStim=goStimBinned, nogoStim=nogoStimBinned, laserStim=laserStimBinned, goStimOnSecs=stimOnSamples[goTrialIndices]/sRate, goStimOffSecs=stimOffSamples[goTrialIndices]/sRate, nogoStimOnSecs=stimOnSamples[nogoTrialIndices]/sRate, nogoStimOffSecs=stimOffSamples[nogoTrialIndices]/sRate, laserStimOnSecs=laserOnsetSamples/sRate, laserStimOffSecs=laserOnsetSamples/sRate+laserDuration)
    save(timeSeries, file=saveFilename)

    browser()
}

processAll()
rm(processAll)
