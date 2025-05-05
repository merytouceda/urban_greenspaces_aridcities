#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=metalsurban
#SBATCH --output=diamond_bacmet.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=125:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load diamond/2.0.9

DB_FILE="/xdisk/egornish/mtoucedasuarez/databases/bacmet/BacMet2_EXP_database.fasta"
READS_DIR="/xdisk/barberan/mtoucedasuarez/urban/clean_reads"
# for loop over samples:
for sample in `awk '{print $1}' /groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/sample_id_list.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # translate the nucleotide fasta reads into protein sequences to be able to run diamond on them

    diamond blastx --db $DB_FILE --query $READS_DIR/${sample}_clean_1.fastq --outfmt 6 --max-target-seqs 1 --evalue 0.000001 --out /xdisk/barberan/mtoucedasuarez/urban/metals/pre-assembly/${sample}_blastp.out --salltitles --quiet --more-sensitive

done
