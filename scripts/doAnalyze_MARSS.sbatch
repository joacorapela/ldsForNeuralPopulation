#!/bin/csh

#SBATCH --job-name=doAnalyze_MARSS
#SBATCH --output=doAnalyz_MARSS_%A_%a.out
#SBATCH --error=doAnalyz_MARSS_%A_%a.err
#SBATCH --array=1-4
#SBATCH --time=00:01:00
#SBATCH --mem=4G

# Print the task id.
echo "Starting task " $SLURM_ARRAY_TASK_ID

set tasksFile='data/tasks.txt'

set stateDim=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $1}' $tasksFile`)
set stateInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $2}' $tasksFile`)
set obsInputMemorySecs=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $3}' $tasksFile`)
set initialCondMethod=(`awk -v myNR=$SLURM_ARRAY_TASK_ID 'NR==myNR {print $4}' $tasksFile`)

echo Rscript $stateDim $stateInputMemorySecs $obsInputMemorySecs $initialCondMethod

echo "Done with task " $SLURM_ARRAY_TASK_ID

