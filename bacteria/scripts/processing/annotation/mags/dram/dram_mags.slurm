#!/bin/bash -l
#SBATCH --job-name=dram-mags
#SBATCH --account=gornish
#SBATCH --output=output_dram.txt
#SBATCH --error=errors_dram.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate DRAM


DRAM.py annotate -i "/xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins/bins_for_virhostmatcher/*.fa" -o /xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins/bins_for_virhostmatcher/annotation
DRAM.py distill -i /xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins/bins_for_virhostmatcher/annotation/annotations.tsv -o /xdisk/barberan/mtoucedasuarez/landuse/mags/bin_refinement/all_bins/bins_for_virhostmatcher/annotation/genome_summaries
