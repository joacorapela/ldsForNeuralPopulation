#!/bin/csh

setenv modelSelectionFile 'data/v1Shaft1_modelSelection.txt'
setenv configFilename 'data/v1Shaft1_estimation.ini'
setenv modelsLogFilename 'log/v1Shaft1Models.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=doAnalyze_MARSS \
--output=slurmOutputs/doAnalyze_MARSS_%A_%a.out \
--error=slurmOutputs/doAnalyze_MARSS_%A_%a.err \
--time=02:00:00 \
--mem=4G \
--array=$arrayOpt \
./doAnalyze_MARSS.sbatch 
