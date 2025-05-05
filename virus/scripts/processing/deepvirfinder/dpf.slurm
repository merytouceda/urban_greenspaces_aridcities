#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=deepvirfinder1
#SBATCH --output=deepvirfinder1.out
#SBATCH --account=laurameredith
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=16:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
#module load anaconda
#source ~/.bashrc && conda activate
#conda activate deepvirfinder
#module load tensorflow/nvidia/2.9.1

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
cd /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus
source config.sh

# array of jobs: 
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

singularity exec --nv /contrib/singularity/shared/sarawillis/biotensor.sif python3 $DVF_BIN/dvf.py -i $IN_DIR/$SMPLE -o $OUT_PATH  -l 1500 -c 8

   
