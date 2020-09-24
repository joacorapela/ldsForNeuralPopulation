create_FA_MARSS <- function(observations, stateDim, stateInputs, stateOffsetType, stateCovType, obsInputs, obsOffsetType, obsCovType, initialStateMeanType, initialStateCovType, maxIter, kfFunc="MARSSkfas", silentLevel=2, controlFA=list(trace=TRUE, nstart=5), controlMARSS=list(trace=1, safe=FALSE), silentMARSS=2) {
    # observations \in obsDim x nObs
    obsDim <- nrow(observations)
    nObs <- ncol(observations)

    # begin estimate initial conditions
    dataForFA <- t(as.matrix(observations))
    initialConds <- estimateKFInitialCondFA(z=dataForFA, nFactors=stateDim, control=controlFA)
    # end estimate initial conditions

    # begin create model
    B1List <- c()
    for(j in 1:stateDim) {
        for(i in 1:stateDim) {
            B1List <- c(B1List, list(sprintf("b%02d%02d", i, j)))
        }
    }
    B1  <- matrix(B1List, nrow=stateDim)
    U1  <- stateOffsetType
    Q1  <- stateCovType
    Z1List <- c()
    for(j in 1:stateDim) {
        for(i in 1:obsDim) {
            Z1List <- c(Z1List, list(sprintf("z%02d%02d", i, j)))
        }
    }
    Z1  <- matrix(Z1List, nrow=obsDim)
    if(!is.na(stateInputs[1])) {
        C1 <- "unconstrained"
        c1 <- stateInputs
    } else {
        C1 <- "zero"
        c1 <- "zero"
    }
    A1  <- obsOffsetType
    R1  <- obsCovType
    pi1 <- initialStateMeanType
    V01 <- initialStateCovType
    if(!is.na(obsInputs[1])) {
        D1 <- "unconstrained"
        d1 <- obsInputs
    } else {
        D1 <- "zero"
        d1 <- "zero"
    }

    model.list <- list(B=B1, U=U1, C=C1, c=c1, Q=Q1, Z=Z1, A=A1, D=D1, d=d1, R=R1, x0=pi1, V0=V01)
    # end create model

    # begin set initial conditions
    B0 <- matrix(as.vector(initialConds$B), ncol=1)
    Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
    R0 <- matrix(initialConds$RDiag, ncol=1)
    controlMARSS <- c(controlMARSS, list(maxit=maxIter))

    inits <- list(B=B0, Z=Z0, R=R0)
    # end set initial conditions

    marsMLE <- MARSS(as.matrix(observations), model=model.list, inits=inits, fun.kf=kfFunc, control=controlMARSS, silent=silentMARSS, fit=FALSE)
    return(marsMLE)
}
