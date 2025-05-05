#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=prodigal_cluster
#SBATCH --output=prodigal_cluster.out
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=168gb
#SBATCH --time=10:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
IN_DIR=/xdisk/barberan/mtoucedasuarez/urban/virus/cluster
OUT_DIR=/xdisk/barberan/mtoucedasuarez/urban/virus/cluster

~/software/Prodigal-2.6.1/prodigal -i $IN_DIR/viral_contigs_rep_seq.fasta -a $OUT_DIR/viral_contigs_cluster_pro.fa -d $OUT_DIR/viral_contigs_cluster_nuc.fa -f gff -p meta -q -m
