#!/bin/csh

#SBATCH --job-name=doAnalyze_MARSS
#SBATCH --output=slurmOutputs/doAnalyze_MARSS_%A_%a.out
#SBATCH --error=slurmOutputs/doAnalyze_MARSS_%A_%a.err
#SBATCH --array=1-4
#SBATCH --time=02:00:00
#SBATCH --mem=8G
# #SBATCH --nodelist=enc1-node13,enc1-node14,enc1-node15,enc2-node7
#SBATCH --exclude=enc3-node1,enc3-node2,enc3-node3,enc3-node4,enc3-node6,enc3-node7

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

set modelSelectionFile='data/v1Shaft1_modelSelection.txt'
set configFilename='data/v1Shaft1_estimation.ini'
set modelsLogFilename='log/v1Shaft1Models.csv'

set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $modelSelectionFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $modelSelectionFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $modelSelectionFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $modelSelectionFile`)

Rscript doAnalyze_MARSS_batch.R --stateDim=$stateDim --stateInputMemorySecs=$stateInputMemorySecs --obsInputMemorySecs=$obsInputMemorySecs --initialCondMethod=$initialCondMethod $configFilename $modelsLogFilename

echo "Done with task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

