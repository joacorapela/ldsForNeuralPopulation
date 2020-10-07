#!/bin/csh

setenv modelSelectionFile 'data/v1Shaft2_modelSelection.txt'
setenv configFilename 'data/v1Shaft2_estimation.ini'
setenv modelsLogFilename 'log/v1Shaft2Models.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=doAnalyze_MARSS \
--output=slurmOutputs/doAnalyze_MARSS_%A_%a.out \
--error=slurmOutputs/doAnalyze_MARSS_%A_%a.err \
--time=02:00:00 \
--mem=4G \
--array=$arrayOpt \
./doAnalyze_MARSS.sbatch 
