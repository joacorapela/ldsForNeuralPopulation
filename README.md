# Code for fitting linear dynamical models with inputs to neural population recordings

## Data preprocessing

To begin data preproscessing you should be in the repository root directory.

1. Create directories and files for the analysis of a new mouse by running `./createMouseDirsAndFiles.csh <mouseName>`, where `<mouseName>` stands for the name of a mouse (e.g., MPV17).

2. Copy the Matlab mouse data file to the directory `data/<mouseName>`. I assume the mouse data filename is `<dataFilenameInfix>.mat`. We will use `<dataFilenameInfix>` in step 6 below.

3. `cd code/scripts`.

4. `matlab`.

5. In the Matlab console type `saveDataSubset(<mouseName>, <dataFilenameInfix>)` to save a subset of the original data in Matlab version 6 format.

6. Edit `../../data/<mouse name>/binLDStimeSeries.ini` and in the caterogry `[filenames]` replace the strings `###` with the `<dataFilenameInfix>` from step 2 above.

7. `Rscript doSaveBinnedLDSTimeSeries.R ../../data/<mouseName>/binLDStimeSeries.ini`. This script will bin spike times.

## Estimation of one model

To begin the estimation of one model you should be in directory `code/scripts'.

8. Create a file with estimation parameters (e.g., `../../data/<mouseName>/v1Shaft1_estimation_DSSSM.ini`)

9. `Rscript doAnalyze_DSSSM_batch.R --stateDim=stateDim --stateInputMemorySecs=stateInputMemorySecs --obsInputMemorySecs=obsInputMemorySecs --initialCondMethod=initialCondMethod binConfigFilename estConfigFilename logFilename`

where initialCondMethod could be FA (factor analysis) or PPCA (probabilistic principal components analysis), binConfigFilename is the file updated in step 6, estConfigFilename is the file created in step 8 and logFilename is the name of an ASCII file where summary information about the estimation of the model will be appended. A line in this file will be of the format:

'<model estimation number> <start time> <train duration> <validation duration> <state dimension> <state memory> <obs memory> <initial conditions method> <log likelihood> <AIC> <cross-validated log likelihood> <estimation elapsed time>'
 (this information includes the random number associated to the estimated model, as well as the log-likelihood and Akaike Information Criterion corresponding to the estimated model).

   For example: `Rscript doAnalyze_DSSSM_batch.R --stateDim=9 --stateInputMemorySecs=0.0 --obsInputMemorySecs=0.6 --initialCondMethod=FA ../../data/<mouseName>/binLDStimeSeries.ini ../../data/<mouseName>/v1Shaft1_estimation_DSSSM.ini ../../log/<mouseName>/v1Shaft1Models_DSSSM.csv`

   This script will print on the screen something like:

   [1] "1045650, 180.000000, 180.000000, 60.000000, 9, 0.000000, 0.800000, PPCA, -12016.202005, 25842.404010, -4667.808143, 65.294000\n"

indicating that the random number associated with the estimated model was 1045650, the models was estimated from a subset of spikes starting at time 180 sec, it was trained using 180 sec of neurla activity, it was validated using 60 sec of neural activity, stateDim=9, stateInputMemorySecs=0, obsInputMemorySecs=0.6, initialCondMethod=PPCA, log likelihood=-12016.2, Akaike Information Criterion=25842.404010, cross-validated Akaike Information Criterion=-4667.808143, estimation elapsed time=65.3 secs. 

This information is saved in `../../log/<mouseName>/v1Shaft1Models_DSSSM.csv`

## Plotting model parameters

To begin the plotting of model parameters you should be in directory `code/scripts`.

`Rscript doPlotDSSSMestimates.R mouseName estNumber`

where estNumber is the random number assigned to the model you want to plot (e.g., first column in `../../log/<mouseName>/v1Shaft1Models_DSSSM.csv`). 

This script will generate a large number of static png files and dynamic html files in the directory `../../figures/<mouseName>`

