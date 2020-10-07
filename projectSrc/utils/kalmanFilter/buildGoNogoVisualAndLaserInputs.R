
buildGoNogoVisualAndLaserInputs <- function(goStim, nogoStim, laserStim, inputMemory) {
    buildInputsBlock <- function(stim, inputMemory) {
        inputs <- c()
        N <- length(stim)
        for(i in 0:(inputMemory)) {
            inputs <- rbind(inputs, c(rep(0, times=i), stim[1:(N-i)]))
        }
        return(inputs)
    }

    goStimInputs <- buildInputsBlock(stim=goStim, inputMemory=inputMemory)
    nogoStimInputs <- buildInputsBlock(stim=nogoStim, inputMemory=inputMemory)
    laserInputs <- buildInputsBlock(stim=laserStim, inputMemory=inputMemory)
    goLaserInputs <- goStimInputs*laserInputs
    nogoLaserInputs <- nogoStimInputs*laserInputs
    allInputs <- rbind(goStimInputs, nogoStimInputs, goLaserInputs, nogoLaserInputs)
    return(allInputs)
}

