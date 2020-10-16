
require(plotly)
require(RColorBrewer)

processAll <- function() {
    myEMestResNumber <- 87430650
    MARSSestResNumber <- 70878808
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estMetaDataFilenamePattern <- "results/%08d_estimation.ini"

    myEMestResFilename <- sprintf(estResFilenamePattern, myEMestResNumber)
    myEMestMetaDataFilename <- sprintf(estMetaDataFilenamePattern, myEMestResNumber)
    MARSSestMetaDataFilename <- sprintf(estMetaDataFilenamePattern, MARSSestResNumber)
    MARSSestResFilename <- sprintf(estResFilenamePattern, MARSSestResNumber)
    myEMestRes <- get(load(myEMestResFilename))
    myEMestMetaData <- read.ini(myEMestMetaDataFilename)
    MARSSestRes <- get(load(MARSSestResFilename))
    MARSSestMetaData <- read.ini(MARSSestMetaDataFilename)

    fig <- plot_ly(type='scatter', mode='markers')
    fig <- fig%>%add_trace(x=1:length(myEMestRes$estRes$logLik), y=myEMestRes$estRes$logLik, name="myEM")
    fig <- fig%>%add_trace(x=1:length(MARSSestRes$kem$iter.record$logLik), y=MARSSestRes$kem$iter.record$logLik, name="MARSS")
    fig <- fig %>% layout(xaxis=list(title="Iteration Number"), yaxis=list(title="Log Likelihood"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_logLikeVsIterNo.%s", myEMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_logLikeVsIterNo.%s", myEMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(fig)

    fig <- plot_ly(type='bar',
                   x=c("MyEM", "MARSS"),
                   y=c(as.numeric(myEMestMetaData$estimation_summary$elapsedTime), 
                       as.numeric(MARSSestMetaData$estimation_summary$elapsedTime)))
    fig <- fig %>% layout(yaxis=list(title="Elapsed Time"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_elapsedTime.%s", myEMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_elapsedTime.%s", myEMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(fig)

    browser()
}

processAll()

