#!/bin/csh

#SBATCH --job-name=doAnalyze_MARSS
#SBATCH --output=slurmOutputs/doAnalyze_MARSS_%A_%a.out
#SBATCH --error=slurmOutputs/doAnalyze_MARSS_%A_%a.err
#SBATCH --array=1-384
#SBATCH --time=02:00:00
#SBATCH --mem=4G

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

set modelSelectionFile='data/lmShaft2_modelSelection.txt'
set configFilename='data/lmShaft2_estimation.ini'
set modelsLogFilename='log/lmShaft2Models.csv'

set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $modelSelectionFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $modelSelectionFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $modelSelectionFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $modelSelectionFile`)

Rscript doAnalyze_MARSS_batch.R --stateDim=$stateDim --stateInputMemorySecs=$stateInputMemorySecs --obsInputMemorySecs=$obsInputMemorySecs --initialCondMethod=$initialCondMethod $configFilename $modelsLogFilename

echo "Done with task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

