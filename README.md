# Code for fitting linear dynamical models with inputs to neural population recordings

## Installation

1. `git clone https://github.com/joacorapela/ldsForNeuralPopulation`

2. `cd ldsForNeuralPopulation`

3. `conda env create -f environment.yml`

4. `conda activate envLDSneuralPopulation`

5. `./installPackages.csh`

## Data preprocessing

1. Save Matlab data in version 6 format with specific a naming convention. See data/doSaveMatlabDataSubset.m (e.g., save this file in
   scripts/results/matlabData_V6.mat)

2. `cd scripts`

3. Bin time series

    a. In data/binLDStimeSeries.ini configure the sample rate (sRate), bin size
(binSizeSecs), laser duration (laserDuration), Matlab data version 6 filename
(generated in step 1; e.g., matlabDataFilename=results/matlabData_V6.mat) and
filename that will contain the binned data (saveFilename),

    b. `Rscript doSaveBinnedLDSTimeSeries.R data/binLDStimeSeries.ini`

## Estimation of one model

1. Create a file with estimation parameters (e.g., data/v1Shaft1_estimation.ini)

2. `Rscript doAnalyze_MARSS_batch.R --stateDim=stateDim --stateInputMemorySecs=stateInputMemorySecs --obsInputMemorySecs=obsInputMemorySecs --initialCondMethod=initialCondMethod estimationParamsFilename logFilename`

    where initialCondMethod could be FA (factor analysis) or PPCA (probabilistic principal components analysis), estimationParamsFilename is the file created in step 1 and logFilename is the name of an ASCII file where summary information about the estimation of the model will be appended (this information includes the random number associated to the estimated model, as well as the log-likelihood and Akaike Information Criterion corresponding to the estimated model).

   For example: `Rscript doAnalyze_MARSS_batch.R --stateDim=9 --stateInputMemorySecs=0.0 --obsInputMemorySecs=0.6 --initialCondMethod=FA data/v1Shaft1_estimation.ini log/v1Shaft1.log`

   This script will print on the screen something like:

   [1] "70190160, 9, 0.000000, 0.600000, FA, -6564.091021, 15384.182043, 15442.239165, 633.024000\n"

    indicating that the random number associated with the estimated model was 70190160, stateDim=9, stateInputMemorySecs=0, obsInputMemorySecs=0.6, initialCondMethod=FA, log likelihood=-6564.091021, Akaike Information Criterion=15384.182043, Akaike Information Criterion corrected=15442.239165, elapsed time=633.024000 secs. 

    This information is saved in log/v1Shaft1.log

## Estimation of multiple models (for model selection)

1. Create a file with estimation parameters (e.g., data/v1Shaft1_estimation.ini)

2. Create a file describing each model to estimate (e.g., data/viShaft1_modelsSelection.txt)

3. In a Unix script file (e.g., doAnalyze_MARSS_v1Shaft1.csh) set the variables configFilename and modelSelectionFile to the filenames created in 1 and 2, and the variable modelsLogFilename to the name of a file where summary information about each estimated model will be appended

4. Run the Unix script file (e.g., `./doAnalyze_MARSS_v1Shaft1.csh`)

The Unix script file will submit in parallel as many jobs to the cluster as models specified in step 2. The best model can be selected as that which maximises the Akaike information criterion reported in modelsLogFilename (step 4).

## Plotting model parameters

1. Configure orca:

    a. `sudo apt-get install xvfb`

    b. type `which orca`

    c. edit the file resulting from step 2 (e.g., vi /nfs/ghome/live/rapela/anaconda3/envs/r_env3/bin/orca) and replace the line:

        exec /nfs/ghome/live/rapela/anaconda3/envs/r_env3/lib/orca_app/orca --no-sandbox "$@"

        with:

        xvfb-run -a  /nfs/ghome/live/rapela/anaconda3/envs/r_env3/lib/orca_app/orca "$@"                                                  

    d. type `orca --help` and you should see something like:

       Plotly's image-exporting utilities

       Usage: orca [--version] [--help] <command> [<args>]

       Available commands:
       - graph [or plotly-graph, plotly_graph]
           Generates an image of plotly graph from inputted plotly.js JSON attributes.
           For more info, run `orca graph --help`.
       - serve [or server]
           Boots up a server with one route per available export component
           For more info, run `orca serve --help`.


3. `Rscript doPlotMARSSestimates.R estNumber`

    where estNumber is the random number assigned to the model you want to plot (e.g., first column in log/v1Shaft1.log). 

    This script will generate a large number of static png files and dynamic html files.

