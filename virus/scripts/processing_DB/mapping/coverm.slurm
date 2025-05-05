#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=coverm
#SBATCH --output=%x.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --mem=470gb
#SBATCH --time=5:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate coverm
# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

cd /xdisk/barberan/mtoucedasuarez/urban/virus/mapping

# require >=95% identity and 45 bp overlap, generate gene count table

TMPDIR=. coverm contig -b *.bam \
-t 94 \
-m count \
--min-read-percent-identity 0.95 \
--min-read-aligned-length 45 \
-o viral_species_count_table.txt
