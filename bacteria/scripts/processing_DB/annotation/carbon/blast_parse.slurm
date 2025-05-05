#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=urb_c_blast_parse
#SBATCH --output=urb_c_blast_parse.out
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

WORKING_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/urban/carbon"
IN_DIR="/xdisk/barberan/mtoucedasuarez/urban/carbon"
OUT_CSV="/xdisk/barberan/mtoucedasuarez/urban/carbon/urb_c_count.csv"

Rscript $WORKING_DIR/blast_parse.R $IN_DIR $OUT_CSV
