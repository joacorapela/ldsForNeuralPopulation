getPlotOneStepAheadForecasts <- function(time, obs, ytt1, Wtt1, obsToPlot=NA, inputs=NA, xlab="Time (sec)", ylab="Square Root Firing Rate", goStimColor="green", nogoStimColor="red", laserStimColor="blue", stimOpacity=0.2) {

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

    getStimRecs <- function(inputs, ymin, ymax) {
        goDiff <- diff(inputs[1,])
        goOnIndices <- which(goDiff>0)+1
        goOffIndices <- which(goDiff<0)+1
        goOnIndices <- goOnIndices[1:min(length(goOnIndices), length(goOffIndices))]
        goOffIndices <- goOffIndices[1:min(length(goOnIndices), length(goOffIndices))]

        goRecs <- getRecs(onsets=time[goOnIndices], durations=time[goOffIndices]-time[goOnIndices], color=goStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

        nogoDiff <- diff(inputs[2,])
        nogoOnIndices <- which(nogoDiff>0)+1
        nogoOffIndices <- which(nogoDiff<0)+1
        nogoOnIndices <- nogoOnIndices[1:min(length(nogoOnIndices), length(nogoOffIndices))]
        nogoOffIndices <- nogoOffIndices[1:min(length(nogoOnIndices), length(nogoOffIndices))]

        nogoRecs <- getRecs(onsets=time[nogoOnIndices], durations=time[nogoOffIndices]-time[nogoOnIndices], color=nogoStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

        laserDiff <- diff(inputs[3,])
        laserOnIndices <- which(laserDiff>0)+1
        laserOffIndices <- which(laserDiff<0)+1
        laserOnIndices <- laserOnIndices[1:min(length(laserOnIndices), length(laserOffIndices))]
        laserOffIndices <- laserOffIndices[1:min(length(laserOnIndices), length(laserOffIndices))]

        laserRecs <- getRecs(onsets=time[laserOnIndices], durations=time[laserOffIndices]-time[laserOnIndices], color=laserStimColor, opacity=stimOpacity, ymin=ymin, ymax=ymax)

        stimRecs <- c(goRecs, nogoRecs, laserRecs)

        return(stimRecs)
    }

    P <- nrow(obs)
    N <- ncol(obs)
    if(is.na(obsToPlot)) {
        obsToPlot <- 1:P
    }

    stdOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        stdOneStepAheadForecast[,i] <- sqrt(diag(Wtt1[,,i]))
    }

    fig <- plot_ly(type='scatter', mode="markers")
    cols <- brewer.pal(max(3, P), "Set1")
    ymax <- -Inf
    ymin <- +Inf
    for(i in obsToPlot) {
        rgbValues <- col2rgb(cols[i])
        observation <- obs[i,]
        forecastMean <- ytt1[i,]
        cbUpper <- ytt1[i,]+1.96*stdOneStepAheadForecast[i,]
        cbLower <- ytt1[i,]-1.96*stdOneStepAheadForecast[i,]
        ymin <- min(ymin, min(c(observation, cbLower)))
        ymax <- max(ymax, max(c(observation, cbUpper)))
        # observation
        fig <- fig%>%add_trace(x=time, y=obs[i,], mode="markers", name=sprintf("observation[,%d]", i), marker=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.5)))
        # forecast
        fig <- fig%>%add_trace(x=time, y=ytt1[i,], mode="lines", name=sprintf("forecast[,%d]", i), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash="solid"))
        fig <- fig%>%add_trace(x=time, y=ytt1[i,]+1.96*stdOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("forecast[,%d]", i), showlegend=FALSE)
        fig <- fig%>%add_trace(x=time, y=ytt1[i,]-1.96*stdOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("forecast[,%d]", i), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    if(!is.na(inputs[1])) {
        stimRecs <- getStimRecs(inputs=inputs, ymin=ymin, ymax=ymax)
        fig <- fig%>%layout(shapes=stimRecs, xaxis=list(title=xlab), yaxis=list(title=ylab))
    } else {
        fig <- fig%>%layout(xaxis=list(title=xlab), yaxis=list(title=ylab))
    }

    return(fig)
}
