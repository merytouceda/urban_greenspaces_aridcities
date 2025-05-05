#!/bin/bash -l                                                                                                                                                               
#SBATCH --job-name=blast-virus                                                                                                                                           
#SBATCH --account=gornish
#SBATCH --output=output_blast_virusmags.txt                                                                                                                                               
#SBATCH --error=errors_blast_virusmags.txt                                                                                                                                                
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu                                        
#SBATCH --partition=standard                                                                                                                                                 
#SBATCH --time=24:00:00                                                                                                                                                      
#SBATCH --ntasks=94                                                                                                                                                           
#SBATCH --nodes=1                                                                                                                                                            
#SBATCH --mem=470gb


# load modules
module load blast
# execute commands

MAGS_FA="/xdisk/barberan/mtoucedasuarez/urban/mags/single/single_70_10/all_bins.fasta"
VIRAL_SPECIES_FASTA="/xdisk/barberan/mtoucedasuarez/urban/virus/cluster/viral_contigs_rep_seq.fasta"
OUTPUT_FILE="/xdisk/barberan/mtoucedasuarez/urban/virus/host-prediction/blast_mags_virus/blast_virus_mags.out"

makeblastdb -in $MAGS_FA -title "Urban mags" -dbtype nucl -out urban_mags_db 
blastn -db urban_mags_db -query $VIRAL_SPECIES_FASTA -out $OUTPUT_FILE -outfmt 6 -evalue 0.0000000001

