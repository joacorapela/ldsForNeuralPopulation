#!/bin/csh

cd ..
set mouseName = exampleMouse
set region = v1
set shaftNro = 1

set time = 02:00:00
set mem = 4G

setenv modelSelectionFile ../../data/{$mouseName}/{$region}Shaft{$shaftNro}_modelSelection.txt
setenv configFilename ../../data/{$mouseName}/{$region}Shaft{$shaftNro}_estimation_DSSSM.ini
setenv modelsLogFilename ../../log/{$mouseName}/{$region}Shaft{$shaftNro}Models_DSSSM.csv

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=dsSSM \
--output=../../slurmOutputs/{$mouseName}/doAnalyze_DSSSM_%A_%a.out \
--error=../../slurmOutputs/{$mouseName}/doAnalyze_DSSSM_%A_%a.err \
--time=$time \
--mem=$mem \
--array=$arrayOpt \
./doAnalyze_DSSSM.sbatch 
cd -

