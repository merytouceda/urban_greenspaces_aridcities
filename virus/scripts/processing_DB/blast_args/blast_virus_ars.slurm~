#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=blastresfinder
#SBATCH --output=blastresfinder.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --mem=48gb
#SBATCH --time=24:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load blast/2.10.0


blastn -db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/resfinder_db/all.fsa -query /xdisk/barberan/mtoucedasuarez/landuse/virus/cluster/viral_contigs_rep_seq.fasta -outfmt 6 -out /xdisk/barberan/mtoucedasuarez/landuse/virus/blast/viral_species_args_blast.out
