#!/bin/bash -l
#SBATCH --job-name=dramv
#SBATCH --account=gornish
#SBATCH --output=output_dramv.txt
#SBATCH --error=errors_dramv.txt
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

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

OUT_PATH="/xdisk/barberan/mtoucedasuarez/landuse/virus/dramv"
VS2_OUTDIR="/xdisk/barberan/mtoucedasuarez/landuse/virus/dramv/vs2_pre_dramV_viralcatalog.out/for-dramv"

DRAM-v.py annotate -i $VS2_OUTDIR/final-viral-combined-for-dramv.fa -v $VS2_OUTDIR/viral-affi-contigs-for-dramv.tab -o annotation
DRAM-v.py distill -i annotation/annotations.tsv -o annotation/distilled