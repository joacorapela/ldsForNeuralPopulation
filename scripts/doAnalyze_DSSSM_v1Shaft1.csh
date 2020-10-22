#!/bin/csh

setenv modelSelectionFile 'data/v1Shaft1_modelSelection_start360_dur180.txt'
setenv configFilename 'data/v1Shaft1_estimation_DSSSM_start360_dur180.ini'
setenv modelsLogFilename 'log/v1Shaft1Models_DSSSM_start360_dur180.csv'

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=doAnalyze_DSSSM \
--output=slurmOutputs/doAnalyze_DSSSM_%A_%a.out \
--error=slurmOutputs/doAnalyze_DSSSM_%A_%a.err \
--time=02:00:00 \
--mem=4G \
--array=$arrayOpt \
./doAnalyze_DSSSM.sbatch 
