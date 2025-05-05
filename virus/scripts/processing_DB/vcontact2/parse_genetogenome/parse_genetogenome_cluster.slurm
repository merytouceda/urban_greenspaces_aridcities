#!/bin/bash -l
#SBATCH --job-name=g2g_cluster
#SBATCH --account=gornish
#SBATCH --output=output_g2g_cluster.txt
#SBATCH --error=errors_g2g_cluster.txt
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu

module load R

SUBMIT_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/vcontact2/parse_genetogenome"
PRODIGAL_DIR="/xdisk/barberan/mtoucedasuarez/urban/virus/cluster"

Rscript $SLURM_SUBMIT_DIR/parse_genetogenome2.R  $PRODIGAL_DIR/viral_contigs_cluster_nuc.fa viral_contigs_cluster
