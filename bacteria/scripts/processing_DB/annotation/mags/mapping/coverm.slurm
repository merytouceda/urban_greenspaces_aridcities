#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=coverm
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --mem=470gb
#SBATCH --time=24:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate coverm
# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

cd /xdisk/barberan/mtoucedasuarez/urban/mags/mapping

MAGS_DIR="/xdisk/barberan/mtoucedasuarez/urban/mags/single/single_70_10"


# require >=95% identity and 45 bp overlap, generate gene count table
TMPDIR=. coverm genome --bam-files *.bam \
-t 94  \
-o mags_count_table.txt \
--genome-fasta-directory $MAGS_DIR \
--genome-fasta-extension fa
