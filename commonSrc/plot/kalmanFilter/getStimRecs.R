
getStimRecs <- function(time, inputs, inputMemorySamples, ymin, ymax, goStimColor="green", nogoStimColor="red", laserStimColor="blue", stimOpacity=0.2) {

    getRecs <- function(onsets, durations, color, opacity, ymin, ymax) {
        recs <- vector(mode="list", length=length(onsets))
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

    goDiff <- diff(inputs[0*(inputMemorySamples+1)+1,])
    goOnIndices <- which(goDiff>0)+1
    goOffIndices <- which(goDiff<0)+1
    goOnIndices <- goOnIndices[1:min(length(goOnIndices), length(goOffIndices))]
    goOffIndices <- goOffIndices[1:min(length(goOnIndices), length(goOffIndices))]

    goRecs <- getRecs(onsets=time[goOnIndices], durations=time[goOffIndices]-time[goOnIndices], color=goStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

    nogoDiff <- diff(inputs[1*(inputMemorySamples+1)+1,])
    nogoOnIndices <- which(nogoDiff>0)+1
    nogoOffIndices <- which(nogoDiff<0)+1
    nogoOnIndices <- nogoOnIndices[1:min(length(nogoOnIndices), length(nogoOffIndices))]
    nogoOffIndices <- nogoOffIndices[1:min(length(nogoOnIndices), length(nogoOffIndices))]

    nogoRecs <- getRecs(onsets=time[nogoOnIndices], durations=time[nogoOffIndices]-time[nogoOnIndices], color=nogoStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

    laserDiff <- diff(inputs[2*(inputMemorySamples+1)+1,])
    laserOnIndices <- which(laserDiff>0)+1
    laserOffIndices <- which(laserDiff<0)+1
    laserOnIndices <- laserOnIndices[1:min(length(laserOnIndices), length(laserOffIndices))]
    laserOffIndices <- laserOffIndices[1:min(length(laserOnIndices), length(laserOffIndices))]

    laserRecs <- getRecs(onsets=time[laserOnIndices], durations=time[laserOffIndices]-time[laserOnIndices], color=laserStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

    stimRecs <- c(goRecs, nogoRecs, laserRecs)

    return(stimRecs)
}

