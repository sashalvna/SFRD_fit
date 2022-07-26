#!/bin/bash
#SBATCH --job-name=CI7          #job name
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks=1            # Number of cores
#SBATCH --output=/n/home04/lvanson/SFRD_fit/src/data//slurm_out/CI7.out            # output storage file
#SBATCH --error=/n/home04/lvanson/SFRD_fit/src/data//slurm_out/CI7.err             # error storage file
#SBATCH --time=1:00:00              # Runtime in minutes
#SBATCH --mem=120000               # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -p demink,conroy,hernquist,shared
#SBATCH --mail-user=aac.van.son@gmail.com         # Send email to user
#SBATCH --mail-type=FAIL       #
#
#Print some stuff on screen
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
#
#Set variables
export QT_QPA_PLATFORM=offscreen # To avoid the X Display error
#
cd /n/home04/lvanson/SFRD_fit/src/scripts//CosmicIntegration/
#
# Run your job
python test.py 
    