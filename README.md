# Code for the analysis of neural population recordings with linear dynamical systems with inpunts

## Installation

## Data preprocessing

1. Save Matlab data in version 6 format using scripts/doSaveDataSubset.m (e.g.,
   save this file in results/matlabData_V6.mat)

2. Bin time series

    a. In data/binLDStimeSeries.ini configure the sample rate (sRate), bin size
(binSizeSecs), laser duration (laserDuration), Matlab data version 6 filename
(generated in the previous step; e.g.,
matlabDataFilename=results/matlabData_V6.mat) and filename that will contain
the binned data (saveFilename),

    b. Rscript doSaveBinnedLDSTimeSeries.R data/binLDStimeSeries.ini

## Estimate one model

1. Create a file with estimation parameters (e.g., data/v1Shaft1_estimation.ini)

2. Rscript doAnalyze_MARSS_batch.R --stateDim=stateDim --stateInputMemorySecs=stateInputMemorySecs --obsInputMemorySecs=obsInputMemorySecs --initialCondMethod=initialCondMethod estimationParamsFilename logFilename

    where initialCondMethod could be FA (factor analysis) or PPCA (probabilisitc principal components analysis), estimationParamsFilename is the file created in step 1 and logFilename is the name of an ASCII file where summary information about the estimation of the model will be appended (this information includes the random number associated to the estimated model, as well as the log-likelihood and Akaike Information Criterion corresponding to the estimated model).

## Estimate multiple models (for model selection)

1. Create a file with estimation parameters (e.g., data/v1Shaft1_estimation.ini)

2. Create a file describing each model to estimate (e.g., data/viShaft1_modelsSelection.txt)

3. In a Unix script file like doAnalyze_MARSS_v1Shaft1.csh set the variables modelSelection?File and configFilename to the filenames created in 1 and 2

4. In the Unix script file set the variable modelsLogFilename to the name of a file where summary information about each estimated model will be appended

5. Run the Unix script file (e.g., ./doAnalyze_MARSS_v1Shaft1.csh)

The Unix script file will submit in parallel as many jobs to the cluster as models specified in step 2. The best model can be selected as that which maximizes the Akaike information criterion reported in modelsLogFilename (step 4).

## Plot estimated model parameters

1. Rscript doPlotMARSSestimates.R estNumber

    where estNumber is the random number assigne to the model you want to plot. This script will generate a large number of static png files and dynamic html files.

