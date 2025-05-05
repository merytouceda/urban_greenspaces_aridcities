#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_index
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=168gb
#SBATCH --time=24:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

IN_DIR="/xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins"
OUT_DIR="/xdisk/barberan/mtoucedasuarez/landuse/mags/mapping"

# generate index
~/software/bwa/bwa index \
$IN_DIR/all_bins_concatenated.fa \
-p $OUT_DIR/mags_index