#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=blast_mags_orfs-hmrgs
#SBATCH --output=outputr_blast_mags_orfs-C_hmrgs.txt
#SBATCH --error=errors_blast_mags_orfs-C_hmrgs.txt
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
BACMET_FA="/xdisk/egornish/mtoucedasuarez/databases/bacmet/BacMet2_EXP_database.fasta"


# blast dbs
#makeblastdb -in $C_FA -title $C_FA -dbtype prot
#makeblastdb -in $BACMET_FA -title $BACMET_FA -dbtype prot


# for loop over samples:
cd $MAGS_ORFS
for file in *.fa_pro.fa 
do
    
        blastp -db $BACMET_FA -query $file -outfmt 6 -out $OUT_DIR/${file}_hmrgs.out  -max_hsps 1 -evalue 0.0000000001 #-salltitles -quiet -more-sensitive

done
