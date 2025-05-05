#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=blast_mags_orfs-args
#SBATCH --output=outputr_blast_mags_orfs-args_N.txt
#SBATCH --error=errors_blast_mags_orfs-args_N.txt
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=40:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load diamond/2.0.9
module load blast/2.13.0 

MAGS_ORFS="/xdisk/barberan/mtoucedasuarez/urban/mags/orfs_mags"
OUT_DIR="/xdisk/barberan/mtoucedasuarez/urban/mags/blast_dbs_to_mags_orfs"
ARGS_BLAST_DB="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/CARD/protein_fasta_protein_homolog_model.fasta"

#diamond makedb --in $C_FA --db $C_DB
#diamond makedb --in $ARGS_FA --db $ARGS_DB
#diamond makedb --in $METALS_FA --db $METALS_DB

# blast dbs
#makeblastdb -in $ARGS_FA -title $ARGS_BLAST_DB -dbtype prot
#makeblastdb -in $N_FA -title $N_BLAST_DB -dbtype prot

# for loop over samples:
cd $MAGS_ORFS
for file in *.fa_pro.fa 
do
    blastp -db $ARGS_BLAST_DB -query $file -outfmt 6 -out $OUT_DIR/${file}_args.out -max_hsps 1 -evalue 0.0000000001 #-salltitles -quiet -more-sensitive
done
