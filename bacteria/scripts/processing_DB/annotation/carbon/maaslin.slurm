#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=urb_c_maaslin
#SBATCH --output=urb_c_maaslin.out
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
COUNT_TABLE="/xdisk/barberan/mtoucedasuarez/urban/carbon/urb_c_count.csv"
cd /xdisk/barberan/mtoucedasuarez/urban/carbon

Rscript /groups/egornish/mtoucedasuarez/scripts/jobs/urban/carbon/maaslin.R $METADATA $COUNT_TABLE
