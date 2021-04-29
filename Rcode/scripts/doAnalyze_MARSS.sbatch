#!/bin/csh

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $modelSelectionFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $modelSelectionFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $modelSelectionFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $modelSelectionFile`)

Rscript doAnalyze_MARSS_batch.R --stateDim=$stateDim --stateInputMemorySecs=$stateInputMemorySecs --obsInputMemorySecs=$obsInputMemorySecs --initialCondMethod=$initialCondMethod $configFilename $modelsLogFilename

echo "Done with task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

