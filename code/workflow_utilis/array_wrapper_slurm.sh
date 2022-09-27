#!/bin/bash

# Job set up:
#SBATCH --array=1-3
#SBATCH --partition=mig
#SBATCH --nodes=1
#SBATCH --job-name="motifbreak_array"
#SBATCH --account="punim0586"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=61440
#SBATCH --time=5-10:0:00
#SBATCH --mail-user=dvespasiani@student.unimelb.edu.au
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# Load modules:
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load r/4.0.0  

## run job (commands must always be at the end)
## specify directories
wd="/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG"
script_dir="$wd/scripts"

command=`head -n $SLURM_ARRAY_TASK_ID ${script_dir}/motifbreak_commands_job_array.sh | tail -n 1`
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo $command
eval $command
