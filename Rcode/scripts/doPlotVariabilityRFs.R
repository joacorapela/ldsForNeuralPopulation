
require(plotly)
require(RColorBrewer)

getPlotMatrixRows <- function(x, aMatrix, legends, xlab="Delay (sec)", ylab="Value", bgColor=rgb(0.35,0.35,0.35)) {
    fig <- plot_ly(type='scatter', mode="lines+markers")
    cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(nrow(aMatrix))
    for(i in 1:nrow(aMatrix)) {
        fig <- fig%>%add_trace(x=x, y=aMatrix[i,], line=list(color=cols[i]), marker=list(color=cols[i]), mode="lines+markers", name=legends[i])
    }
    fig <- fig%>%layout(xaxis=list(title=xlab), yaxis=list(title=ylab),plot_bgcolor=bgColor)
    return(fig)
}

processAll <- function() {
    DEBUG <- TRUE
    if(!DEBUG) {
        option_list <- list(
            make_option(c("-b", "--bestModelsFilenamePattern"), type="character", default="../../results/%s/%sShaft%s_DSSSM_AIC_bestModelsByStartTime.csv", help="Best models filename pattern"),
            make_option(c("-e", "--estResFilenamePattern"), type="character", default="../../results/%08d_estimation.RData", help="Estimation result filename pattern"),
            make_option(c("-f", "--variabilityFigFilenamePattern"), type="character", default="../../figures/%s/%sShaft%s_DSSSM_%s_neuron%d.%s", help="Figure filename pattern"),
            make_option(c("-s", "--sRate"), type="int", default=10, help="Sample rate"),
            make_option(c("-s", "--stateInputMemoryToPlotSamples"), type="int", default=6, help="State input memory to plot (in samples)"),
        )
        parser <- OptionParser(usage = "%prog [options] mouseName region shaft", option_list=option_list)
        parseRes <- parse_args(parser, positional_arguments=3)
        options <- parseRes$options
        arguments <- parseRes$args

        mouseName <- arguments[1]
        region <- arguments[2]
        shaft <- arguments[3]
        bestModelsFilenamePattern <- options$bestModelsFilenamePattern
        estResFilenamePattern <- options$estResFilenamePattern
        variabilityFigFilenamePattern <- options$variabilityFigFilenamePattern
        sRate <- options$sRate
        stateInputMemoryToPlotSamples <- options$stateInputMemoryToPlotSamples
    } else {
        mouseName <- "VL61"
        region <- "v1"
        shaft <- "1"
        bestModelsFilenamePattern <- "../../results/%s/%sShaft%s_DSSSM_AIC_bestModelsByStartTime.csv"
        estResFilenamePattern <- "../../results/%s/%08d_estimation.RData"
        variabilityFigFilenamePattern <- "../../figures/%s/%sShaft%s_DSSSM_%s_neuron%d.%s"
        sRate <- 10
        stateInputMemoryToPlotSamples <- 6
    }
    bestModelsFilename <- sprintf(bestModelsFilenamePattern, mouseName, region, shaft)
    bestModels <- read.table(file=bestModelsFilename, header=TRUE)
    sortRes <- sort(bestModels$analysisStartTimeSecs, index.return=TRUE)
    sortedBestModels <- bestModels[sortRes$ix,]
    obsInputMemorySamples <- sortedBestModels[1, "obsInputMemorySecs"]*sRate
    stateInputMemorySamples <- sortedBestModels[1, "stateInputMemorySecs"]*sRate
    if(stateInputMemorySamples!=0) {
        stop("Currently only stateInputMemorySamples=0 is supported")
    }
    xD <- (0:obsInputMemorySamples)/sRate
    xC <- (0:stateInputMemoryToPlotSamples)/sRate
    DblockSize <- 1+obsInputMemorySamples
    CblockSize <- 1+stateInputMemorySamples
    legends <- rep(NA, times=nrow(sortedBestModels))
    legendPattern <- "start %d sec"

    estNumber <- sortedBestModels[1, "estNumber"]
    estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
    estRes <- get(load(estResFilename))
    dsSSM <- estRes$dsSSM
    nNeurons <- nrow(dsSSM$Z)

    for(n in 1:nNeurons) {

if(FALSE) {
        # start DVisualStimGo
        Doffset <- 0
        figDesc <- "DVisualStimGo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xD))
            }
            rfs[i,] <- dsSSM$D[n,Doffset+(1:DblockSize)]
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xD, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end DVisualStimGo

        # start DVisualStimNogo
        Doffset <- Doffset+DblockSize
        figDesc <- "DVisualStimNogo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xD))
            }
            rfs[i,] <- dsSSM$D[n,Doffset+(1:DblockSize)]
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xD, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end DVisualStimNogo

        # start DLaserStimGo
        Doffset <- Doffset+DblockSize
        figDesc <- "DLaserStimGo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xD))
            }
            rfs[i,] <- dsSSM$D[n,Doffset+(1:DblockSize)]
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xD, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end DLaserStimGo

        # start DLaserStimNogo
        Doffset <- Doffset+DblockSize
        figDesc <- "DLaserStimNogo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xD))
            }
            rfs[i,] <- dsSSM$D[n,Doffset+(1:DblockSize)]
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xD, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end DLaserStimNogo
}
        # start ZBCVisualStimGo
        Coffset <- 0
        figDesc <- "ZBCVisualStimGo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xC))
            }
            Bpow <- diag(rep(1, times=nrow(dsSSM$B)))
            for(l in 0:(ncol(rfs)-1)) {
                rfs[i,l+1] <- dsSSM$Z[n,]%*%Bpow%*%dsSSM$C[,Coffset+(1:CblockSize)]
                Bpow <- Bpow%*%dsSSM$B
            }
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xC, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end ZBCVisualStimGo

        # start ZBCVisualStimNogo
        Coffset <- Coffset+CblockSize
        figDesc <- "ZBCVisualStimNogo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xC))
            }
            Bpow <- diag(rep(1, times=nrow(dsSSM$B)))
            for(l in 0:(ncol(rfs)-1)) {
                rfs[i,l+1] <- dsSSM$Z[n,]%*%Bpow%*%dsSSM$C[,Coffset+(1:CblockSize)]
                Bpow <- Bpow%*%dsSSM$B
            }
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xC, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end ZBCVisualStimNogo

        # start ZBCLaserStimGo
        Coffset <- Coffset+CblockSize
        figDesc <- "ZBCLaserStimGo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xC))
            }
            Bpow <- diag(rep(1, times=nrow(dsSSM$B)))
            for(l in 0:(ncol(rfs)-1)) {
                rfs[i,l+1] <- dsSSM$Z[n,]%*%Bpow%*%dsSSM$C[,Coffset+(1:CblockSize)]
                Bpow <- Bpow%*%dsSSM$B
            }
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xC, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end ZBCLaserStimGo

        # start ZBCLaserStimNogo
        Coffset <- Coffset+CblockSize
        figDesc <- "ZBCLaserStimNogo"
        show(sprintf("%s, neuron %d", figDesc, n))
        for(i in 1:nrow(sortedBestModels)) {
            estNumber <- sortedBestModels[i, "estNumber"]
            estResFilename <- sprintf(estResFilenamePattern, mouseName, estNumber)
            estRes <- get(load(estResFilename))
            dsSSM <- estRes$dsSSM
            if(i==1) {
                rfs <- matrix(NA, nrow=nrow(sortedBestModels), ncol=length(xC))
            }
            Bpow <- diag(rep(1, times=nrow(dsSSM$B)))
            for(l in 0:(ncol(rfs)-1)) {
                rfs[i,l+1] <- dsSSM$Z[n,]%*%Bpow%*%dsSSM$C[,Coffset+(1:CblockSize)]
                Bpow <- Bpow%*%dsSSM$B
            }
            legends[i] <- sprintf(legendPattern, sortedBestModels[i, "analysisStartTimeSecs"])
        }
        fig <- getPlotMatrixRows(x=xC, aMatrix=rfs, legends=legends)
        htmlFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "html")
        pngFilename <- sprintf(variabilityFigFilenamePattern, mouseName, region, shaft, figDesc, n, "png")
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
        orca(p=fig, file=pngFilename)
        # print(fig)
        # end ZBCLaserStimNogo
    }
    browser()
}

processAll()
