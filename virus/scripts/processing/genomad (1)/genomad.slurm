#!/bin/bash

#SBATCH --job-name=genomad_clustered
#SBATCH -o %x.out
#SBATCH -e %x.err 
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=470gb
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=barberan

source ~/.bashrc
mamba activate
mamba activate genomad

INPUT_FILE="/xdisk/egornish/mtoucedasuarez/barberan_xdisk_0224/landuse/virus/cluster/viral_contigs_rep_seq.fasta"
OUTPUT_DIR="/xdisk/egornish/mtoucedasuarez/barberan_xdisk_0224/landuse/virus/genomad/clustered"
DATABASE="/xdisk/egornish/mtoucedasuarez/databases/genomad_db"

genomad end-to-end --cleanup $INPUT_FILE $OUTPUT_DIR $DATABASE
