#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=LU_carbon_strict
#SBATCH --output=LU_carbon_strict.out
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
module load diamond/2.0.9
source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/carbon/config.sh

export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
QUERY=$READS_DIR/fasta/${SMPLE}_R1_paired.fa 

# array command
diamond blastx --db $DB_DIR/CAZyDB.09242021.fa --query $QUERY --outfmt 6 --max-target-seqs 1 --evalue 0.000000001 --out $OUT_PATH/${SMPLE}_blastx_strict.out --salltitles --quiet --more-sensitive



