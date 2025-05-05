#!/bin/bash

#SBATCH --job-name=mmseqs2
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=470gb
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=barberan


module load anaconda
source ~/.bashrc && conda activate
conda activate mmseqs2

cd /groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/cluster
# cluster based on 95% identity and 90% overlap, we use --cov-mode 1 to being related to cd-hit
# See https://mmseqs.com/latest/userguide.pdf
mmseqs easy-cluster /xdisk/barberan/mtoucedasuarez/urban/virus/final/fasta/all_viral_contigs.fa viral_contigs tmp --min-seq-id 0.95 --cov-mode 1 -c 0.7 --threads 94
