#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=virhostmatcher
#SBATCH --output=%x.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=168gb
#SBATCH --time=24:00:00

module load python/3.6


python /groups/barberan/software/VirHostMatcher/vhm.py -v /xdisk/barberan/mtoucedasuarez/landuse/virus/final/fasta/ -b /xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins/bins_for_virhostmatcher -o /xdisk/barberan/mtoucedasuarez/landuse/virus/virhostmatcher -t /xdisk/barberan/mtoucedasuarez/landuse/mags/mags_tax.txt 
