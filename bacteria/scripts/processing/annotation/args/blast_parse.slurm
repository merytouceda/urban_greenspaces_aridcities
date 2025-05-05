#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=LU_card_blast_parse
#SBATCH --output=LU_cards_blast_parse.out
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
module load R/4.2.2 

WORKING_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args"
IN_DIR="/xdisk/barberan/mtoucedasuarez/landuse/args/card/pre-assembly"
OUT_CSV="/xdisk/barberan/mtoucedasuarez/landuse/args/card/lu_card_count.csv"

Rscript $WORKING_DIR/blast_parse.R $IN_DIR $OUT_CSV
