#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=ncyc_100_fasta_len
#SBATCH --output=ncyc_100_fasta_len.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=20:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load python/3.9/3.9.10

INPUT_DB="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/nitrogen/NCyc/data/NCyc_100.faa"
OUT_FILE="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/nitrogen/NCyc/data/ncyc_100_len.csv"

python3 fasta_len.py $INPUT_DB $OUT_FILE
