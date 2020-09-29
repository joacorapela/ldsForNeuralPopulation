getPlotSmoothedStates <- function(time, xtT, VtT, stateInputMemorySamples, inputs=NA, statesToPlot=NA, xlab="Time (sec)", ylab="p(x|y_1,...,yN)", goStimColor="green", nogoStimColor="red", laserStimColor="blue", stimOpacity=0.2) {

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
        rgbaColorName <- sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1)
        fig <- fig%>%add_trace(x=time, y=mean, mode="lines+markers", name=sprintf("state%d", i), line=list(color=rgbaColorName, dash="solid"), marker=list(color=rgbaColorName, symbol="asterisk-open"))
        fig <- fig%>%add_trace(x=time, y=cbUpper, mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("state[,%d]", i), showlegend=FALSE)
        fig <- fig%>%add_trace(x=time, y=cbLower, mode="lines", line=list(color="rgba(0,0,0,0)"), name=sprintf("state[,%d]", i), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    fig <- fig%>%layout(xaxis=list(title=xlab), yaxis=list(title=ylab))
    if(!is.na(inputs[1])) {
        stimRecs <- getStimRecs(time=time, inputs=inputs, inputMemorySamples=stateInputMemorySamples, ymin=ymin, ymax=ymax, goStimColor=goStimColor, nogoStimColor=nogoStimColor, laserStimColor=laserStimColor, stimOpacity=stimOpacity)
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
