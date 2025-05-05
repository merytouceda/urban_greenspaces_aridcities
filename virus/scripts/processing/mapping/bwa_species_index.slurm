#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_index_5000
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

IN_DIR="/xdisk/egornish/mtoucedasuarez/barberan_xdisk_0224/landuse/virus/cluster"
OUT_DIR="/xdisk/egornish/mtoucedasuarez/barberan_xdisk_0224/landuse/virus/mapping/5000"

# generate index
~/software/bwa/bwa index \
$IN_DIR/viral_contigs_5000_rep_seq.fasta \
-p $OUT_DIR/species_index_5000
