
require(plotly)
require(ini)
require(RColorBrewer)

processAll <- function() {
    dsSSMestResNumber <- 34633421
    MARSSestResNumber <- 70878808
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estMetaDataFilenamePattern <- "results/%08d_estimation.ini"

    dsSSMestResFilename <- sprintf(estResFilenamePattern, dsSSMestResNumber)
    dsSSMestMetaDataFilename <- sprintf(estMetaDataFilenamePattern, dsSSMestResNumber)
    MARSSestMetaDataFilename <- sprintf(estMetaDataFilenamePattern, MARSSestResNumber)
    MARSSestResFilename <- sprintf(estResFilenamePattern, MARSSestResNumber)
    dsSSMestRes <- get(load(dsSSMestResFilename))
    dsSSMestMetaData <- read.ini(dsSSMestMetaDataFilename)
    MARSSestRes <- get(load(MARSSestResFilename))
    MARSSestMetaData <- read.ini(MARSSestMetaDataFilename)

    fig <- plot_ly(type='scatter', mode='markers')
    fig <- fig%>%add_trace(x=1:length(dsSSMestRes$dsSSM$logLik), y=dsSSMestRes$dsSSM$logLik, name="dsSSM")
    fig <- fig%>%add_trace(x=1:length(MARSSestRes$kem$iter.record$logLik), y=MARSSestRes$kem$iter.record$logLik, name="MARSS")
    fig <- fig %>% layout(xaxis=list(title="Iteration Number"), yaxis=list(title="Log Likelihood"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_dsSSMvsMARSS_logLikeVsIterNo.%s", dsSSMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_dsSSMvsMARSS_logLikeVsIterNo.%s", dsSSMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    fig <- plot_ly(type='bar',
                   x=c("DSSSM", "MARSS"),
                   y=c(as.numeric(dsSSMestMetaData$estimation_summary$elapsedTime), 
                       as.numeric(MARSSestMetaData$estimation_summary$elapsedTime)))
    fig <- fig %>% layout(yaxis=list(title="Elapsed Time"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_dsSSMvsMARSS_elapsedTime.%s", dsSSMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_dsSSMvsMARSS_elapsedTime.%s", dsSSMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    browser()
}

processAll()

