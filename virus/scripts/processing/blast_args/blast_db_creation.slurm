#!/bin/bash -l
#SBATCH --job-name=blast_db
#SBATCH --account=gornish
#SBATCH --output=output_blast_db.txt
#SBATCH --error=errors_blast_db.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu


# We are going to create a blast db with our viral species (output from clustering)

module load blast/2.10.0

IN_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/cluster"

makeblastdb -in $IN_DIR/viral_contigs_rep_seq.fasta -out $IN_DIR/viral_species_blast_db -dbtype nucl -title viral_species_landuse -parse_seqids
