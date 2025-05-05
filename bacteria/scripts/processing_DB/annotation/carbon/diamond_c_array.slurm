#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=urb_c_strict
#SBATCH --output=urb_c_strict.out
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
source /groups/egornish/mtoucedasuarez/scripts/jobs/urban/carbon/config.sh

export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
QUERY=$READS_DIR/${SMPLE}_clean_1.fastq


diamond blastx --db $DB_DIR/CAZyDB.09242021.fa --query $QUERY --outfmt 6 --out $OUT_DIR/${SMPLE}_blastx_strict.out -k 1 --evalue 0.000000001  --salltitles  --quiet --more-sensitive
