#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=urb_card_maaslin
#SBATCH --output=urb_card_maaslin.out
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

METADATA="/xdisk/barberan/mtoucedasuarez/urban/metadata.csv"
COUNT_TABLE="/xdisk/barberan/mtoucedasuarez/urban/args/card/urb_card_count.csv"
cd /xdisk/barberan/mtoucedasuarez/urban/args/card

Rscript /groups/egornish/mtoucedasuarez/scripts/jobs/urban/args/card/maaslin.R $METADATA $COUNT_TABLE
