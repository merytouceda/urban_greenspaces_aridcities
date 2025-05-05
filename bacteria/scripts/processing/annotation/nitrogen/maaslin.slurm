#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=lu_n_maaslin
#SBATCH --output=lu_n_maaslin.out
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

METADATA="/xdisk/barberan/mtoucedasuarez/landuse/LU_metadata.csv"
COUNT_TABLE="/xdisk/barberan/mtoucedasuarez/landuse/nitrogen/lu_n_count.csv"
cd /xdisk/barberan/mtoucedasuarez/landuse/nitrogen

Rscript /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/nitrogen/maaslin.R $METADATA $COUNT_TABLE
