#!/bin/csh

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

set analysisStartTimeSecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $modelSelectionFile`)
set trainDurSecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $modelSelectionFile`)
set validationDurSecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $modelSelectionFile`)
set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $modelSelectionFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $5}' $modelSelectionFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $6}' $modelSelectionFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $7}' $modelSelectionFile`)

Rscript doAnalyze_DSSSM_batch.R --analysisStartTimeSecs=$analysisStartTimeSecs --trainDurSecs=$trainDurSecs --validationDurSecs=$validationDurSecs --stateDim=$stateDim --stateInputMemorySecs=$stateInputMemorySecs --obsInputMemorySecs=$obsInputMemorySecs --initialCondMethod=$initialCondMethod $configFilename $modelsLogFilename

echo "Done with task " $SLURM_ARRAY_TASK_ID " on " $SLURMD_NODENAME

