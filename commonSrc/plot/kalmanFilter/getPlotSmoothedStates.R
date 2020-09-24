getPlotSmoothedStates <- function(time, xtT, VtT, inputs=NA, statesToPlot=NA, xlab="Time (sec)", ylab="p(x|y_1,...,yN)", goStimColor="green", nogoStimColor="red", laserStimColor="blue", stimOpacity=0.2) {
    # inputs[1,]: go
    # inputs[2,]: nogo
    # inputs[3,]: laser

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

    M <- nrow(xtT)
    N <- ncol(xtT)
    if(!is.na(inputs[1])) {
        nInputs <- nrow(inputs)

    } else {
        nInputs <- 0
    }
    if(is.na(statesToPlot[1])) {
        statesToPlot <- 1:M
    }

    stds <- matrix(NA, nrow=M, ncol=N)
    for(i in 1:N) {
        stds[,i] <- sqrt(diag(VtT[,,i]))
    }

    fig <- plot_ly(type='scatter', mode="markers")
    # cols <- brewer.pal(max(3, M+nInputs), "Set1")
    cols <- colorRampPalette(brewer.pal(9, "Set1"))(max(3, M+nInputs))
    ymax <- -Inf
    ymin <- +Inf
    for(i in statesToPlot) {
        rgbValues <- col2rgb(cols[i])
        mean <- xtT[i,]
        std <- stds[i,]
        cbUpper <- mean+1.96*std
        cbLower <- mean-1.96*std
        ymax <- max(ymax, max(cbUpper))
        ymin <- min(ymin, min(cbLower))
        fig <- fig%>%add_trace(x=time, y=mean, mode="lines", name=sprintf("state%d", i), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash="solid"))
        fig <- fig%>%add_trace(x=time, y=cbUpper, mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("forecast[,%d]", i), showlegend=FALSE)
        fig <- fig%>%add_trace(x=time, y=cbLower, mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("forecast[,%d]", i), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    fig <- fig%>%layout(xaxis=list(title=xlab), yaxis=list(title=ylab))
    if(!is.na(inputs[1])) {
        stimRecs <- getStimRecs(inputs=inputs, ymin=ymin, ymax=ymax)
        fig <- fig%>%layout(shapes=stimRecs)
    }
#     if(!is.na(inputs[1])) {
#         for(i in 1:nrow(inputs)) {
#             rgbValues <- col2rgb(cols[M+i])
#             fig <- fig%>%add_trace(x=time, y=inputs[i,], mode="lines", name=sprintf("input%d", i), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash="solid"))
#         }
#     }

    return(fig)
}
