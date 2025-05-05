#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_index
#SBATCH --output=%x.out
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

IN_DIR="/xdisk/barberan/mtoucedasuarez/urban/virus/cluster"
OUT_DIR="/xdisk/barberan/mtoucedasuarez/urban/virus/mapping"

# generate index
~/software/bwa/bwa index \
$IN_DIR/viral_contigs_rep_seq.fasta \
-p $OUT_DIR/species_index
