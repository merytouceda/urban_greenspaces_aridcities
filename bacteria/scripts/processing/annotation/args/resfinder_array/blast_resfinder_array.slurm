#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=blast_resfinder
#SBATCH --output=blast_resfinder_array.out
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
# Load modules here
# --------------------------------------------------
module load blast/2.13.0 
source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/resfinder_array/config.sh

export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
QUERY=$READS_DIR/${SAMPLE}_R1_paired.fa

# array command
blastn -db $DB_DIR/all.fsa -query $QUERY -outfmt 6 -max_target_seqs 1 -evalue 0.000001 -out $OUT_PATH/${SAMPLE}_blast.out

