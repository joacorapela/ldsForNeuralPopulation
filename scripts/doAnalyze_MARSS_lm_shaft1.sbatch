#!/bin/csh

#SBATCH --job-name=doAnalyze_MARSS
#SBATCH --output=slurmOutputs/doAnalyze_MARSS_%A_%a.out
#SBATCH --error=slurmOutputs/doAnalyze_MARSS_%A_%a.err
#SBATCH --array=1-180
#SBATCH --time=02:00:00
#SBATCH --mem=4G

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID

set tasksFile='data/tasks_lm_shaft1.txt'
set configFilename='data/00000003_estimation.ini'
set modelsLogFilename='log/lmShaft1Models.csv'

set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $tasksFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $tasksFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $tasksFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $tasksFile`)

Rscript doAnalyze_MARSS_batch.R --stateDim=$stateDim --stateInputMemorySecs=$stateInputMemorySecs --obsInputMemorySecs=$obsInputMemorySecs --initialCondMethod=$initialCondMethod $configFilename $modelsLogFilename

echo "Done with task " $SLURM_ARRAY_TASK_ID

