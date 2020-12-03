#!/bin/csh

cd ..
set cellName = VL61
set region = lm
set shaftNro = 1

set time = 02:00:00
set mem = 4G

setenv modelSelectionFile ../../data/{$cellName}/{$region}Shaft{$shaftNro}_modelSelection.txt
setenv configFilename ../../data/{$cellName}/{$region}Shaft{$shaftNro}_estimation_DSSSM.ini
setenv modelsLogFilename ../../log/{$cellName}/{$region}Shaft{$shaftNro}Models_DSSSM.csv

set arrayOpt=1-`wc -l < $modelSelectionFile`

sbatch \
--job-name=dsSSM \
--output=../../slurmOutputs/{$cellName}/doAnalyze_DSSSM_%A_%a.out \
--error=../../slurmOutputs/{$cellName}/doAnalyze_DSSSM_%A_%a.err \
--time=$time \
--mem=$mem \
--array=$arrayOpt \
./doAnalyze_DSSSM.sbatch 
cd -
