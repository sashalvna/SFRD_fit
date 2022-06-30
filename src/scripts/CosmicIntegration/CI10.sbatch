#!/bin/bash
#SBATCH --job-name=CI10          #job name
#SBATCH --nodes=1             # Number of nodes
#SBATCH --ntasks=1            # Number of cores
#SBATCH --output=../../data/slurm_out/CI10.out            # output storage file
#SBATCH --error=../../data/slurm_out/CI10.err             # error storage file
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
# Run your job
python FastCosmicIntegration.py  --path /n/home04/lvanson/sfr_copy/SFRD_fit/src/data/ --filename COMPAS_Output_wWeights.h5 --outfname Rate_info.hdf5 --mu0 0.025 --muz -0.049 --sigma0 1.129 --sigmaz 0.048 --alpha -0.9 --aSF 0.017 --bSF 1.481 --cSF 4.452 --dSF 5.913 --weight mixture_weight --zstep 0.01 --sens O3 --m1min 10. --dco_type BBH --redshiftBinSize 0.05 --maxzdet 0.5 > ../../data/slurm_out/CI10.log
