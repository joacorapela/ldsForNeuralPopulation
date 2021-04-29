processAll <- function() {
    DEBUG <- TRUE
    if(!DEBUG) {
        option_list <- list(
            make_option(c("-m", "--modelsLogFilenamePattern"), type="character", default="../../log/%s/%sShaft%dModels_DSSSM.csv", help="Models log filename pattern"),
            make_option(c("-f", "--bestModelsFilenamePattern"), type="character", default="../../log/%s/%sShaft%d_DSSSM_AIC_bestModelsByStartTime.csv", help="Best models filename pattern"),
            make_option(c("-f", "--bestModelsFigFilenamePattern"), type="character", default="../../figures/%s/%sShaft%d_DSSSM_AIC_bestModelsByStartTime.png", help="Best models figure filename pattern"),
        )
        parser <- OptionParser(usage = "%prog [options] mouseName region shaft", option_list=option_list)
        parseRes <- parse_args(parser, positional_arguments=3)
        options <- parseRes$options
        arguments <- parseRes$args

        mouseName <- arguments[1]
        region <- arguments[2]
        shaft <- as.numeric(arguments[3])
        modelsLogFilenamePattern <- options$modelsLogFilenamePattern
        bestModelsFilenamePattern <- options$bestModelsFilenamePattern
        bestModelsFigFilenamePattern <- options$bestModelsFigFilenamePattern
    } else {
        mouseName <- "VL61"
        region <- "v1"
        shaft <- 1
        modelsLogFilenamePattern <- "../../log/%s/%sShaft%dModels_DSSSM.csv"
        bestModelsFilenamePattern <- "../../results/%s/%sShaft%d_DSSSM_AIC_bestModelsByStartTime.csv"
        bestModelsFigFilenamePattern <- "../../figures/%s/%sShaft%d_DSSSM_AIC_bestModelsByStartTime.png"
    }

    modelsLogFilename <- sprintf(modelsLogFilenamePattern, mouseName, region, shaft)
    bestModelsFilename <- sprintf(bestModelsFilenamePattern, mouseName, region, shaft)
    bestModelsFigFilename <- sprintf(bestModelsFigFilenamePattern, mouseName, region, shaft)

    modelsLog <- read.table(modelsLogFilename, sep=",")
    colnames(modelsLog) <- c("estNumber", "analysisStartTimeSecs", "trainDurSecs", "validationDurSecs", "stateDim", "stateInputMemorySecs", "obsInputMemorySecs", "initialCondMethod", "logLike", "AIC", "cvLogLike", "elapsedTime")
    startTimes <- unique(modelsLog$analysisStartTimeSecs)
    bestModels <- c()
    for(i in 1:length(startTimes)) {
        modelsLogForStartTime <- modelsLog[modelsLog$analysisStartTimeSecs==startTimes[i],]
        selectedModelIndex <- which.min(modelsLogForStartTime$AIC)
        bestModels <- rbind(bestModels, modelsLogForStartTime[selectedModelIndex,])
    }
    write.table(bestModels, file=bestModelsFilename, row.names=FALSE)

    res <- sort(bestModels$analysisStartTimeSecs, index.return=TRUE)
    png(bestModelsFigFilename)
    plot(bestModels$analysisStartTimeSecs[res$ix], bestModels$stateDim[res$ix], type="b", xlab="Start Time (sec)", ylab="Optinmal State Dimension")
    dev.off()

    browser()
}

processAll()

rm(processAll)
