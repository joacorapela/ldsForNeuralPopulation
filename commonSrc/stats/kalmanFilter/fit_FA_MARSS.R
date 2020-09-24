fit_FA_MARSS <- function(observations, stateDim, stateInputs, stateOffsetType, stateCovType, obsInputs, obsOffsetType, obsCovType, initialStateMeanType, initialStateCovType, maxIter, kfFunc, silentLevel=2, controlFA=list(trace=TRUE, nstart=5)) {
    MLEobj <- create_FA_MARSS(observations=observations,
                              stateDim=stateDim,
                              stateInputs=stateInputs,
                              stateOffsetType=stateOffsetType,
                              stateCovType=stateCovType,
                              obsInputs=obsInputs,
                              obsOffsetType=obsOffsetType,
                              obsCovType=obsCovType,
                              initialStateMeanType=initialStateMeanType,
                              initialStateCovType=initialStateCovType,
                              maxIter=maxIter,
                              kfFunc=kfFunc,
                              silentLevel=silentLevel,
                              controlFA=controlFA)
    kem <- MARSSkem(MLEobj=MLEobj)
    return(kem)
}
