#!/bin/csh

setenv modelSelectionFile 'data/lmShaft1_modelSelection.txt'
setenv configFilename 'data/lmShaft1_estimation.ini'
setenv modelsLogFilename 'log/lmShaft1Models.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=doAnalyze_MARSS \
--output=slurmOutputs/doAnalyze_MARSS_%A_%a.out \
--error=slurmOutputs/doAnalyze_MARSS_%A_%a.err \
--time=02:00:00 \
--mem=4G \
--array=$arrayOpt \
./doAnalyze_MARSS.sbatch 
