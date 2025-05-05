#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=urb_fargene
#SBATCH --output=urb_fargene.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=16:00:00


# --------------------------------------------------
# Load modules
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate seqkit

# --------------------------------------------------
# Execute
# --------------------------------------------------
cd /groups/egornish/mtoucedasuarez/scripts/jobs/urban/args/fargene 
source config.sh


export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

fargene --infiles $READS_DIR/${SAMPLE}_contigs.faa --hmm-model class_a --no-quality-filtering --no-orf-prediction