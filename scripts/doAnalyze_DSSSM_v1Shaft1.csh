#!/bin/csh

setenv modelSelectionFile 'data/v1Shaft1_modelSelection.txt'
setenv configFilename 'data/v1Shaft1_estimation_DSSSM.ini'
setenv modelsLogFilename 'log/v1Shaft1Models_DSSSM.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=dsSSM \
--output=slurmOutputs/doAnalyze_DSSSM_%A_%a.out \
--error=slurmOutputs/doAnalyze_DSSSM_%A_%a.err \
--time=02:00:00 \
--mem=4G \
--array=$arrayOpt \
./doAnalyze_DSSSM.sbatch 
