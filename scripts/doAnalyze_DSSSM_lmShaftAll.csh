#!/bin/csh

setenv modelSelectionFile 'data/lmShaftAll_modelSelection.txt'
setenv configFilename 'data/lmShaftAll_estimation_DSSSM.ini'
setenv modelsLogFilename 'log/lmShaftAllModels_DSSSM.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=doAnalyze_MARSS \
--output=slurmOutputs/doAnalyze_MARSS_%A_%a.out \
--error=slurmOutputs/doAnalyze_MARSS_%A_%a.err \
--time=4:00:00 \
--mem=6G \
--array=$arrayOpt \
./doAnalyze_DSSSM.sbatch 
