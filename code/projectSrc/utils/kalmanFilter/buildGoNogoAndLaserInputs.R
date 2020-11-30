
buildGoNogoAndLaserInputs <- function(goStim, nogoStim, laserStim, inputMemory) {
    buildInputsBlock <- function(stim, inputMemory) {
        inputs <- c()
        N <- length(stim)
        for(i in 0:(inputMemory)) {
            inputs <- rbind(inputs, c(rep(0, times=i), stim[1:(N-i)]))
        }
        return(inputs)
    }

    goInputs <- buildInputsBlock(stim=goStim, inputMemory=inputMemory)
    nogoInputs <- buildInputsBlock(stim=nogoStim, inputMemory=inputMemory)
    laserInputs <- buildInputsBlock(stim=laserStim, inputMemory=inputMemory)
    allInputs <- rbind(goInputs, nogoInputs, laserInputs)                                                                              
    return(allInputs)
}   

