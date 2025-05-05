#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=cardLU
#SBATCH --output=diamondacard.out
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
#source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/config.sh

#export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
#$QUERY=$IN_DIR/$SMPLE
# array command
#diamond blastp --db $DB_DIR/protein_fasta_protein_homolog_model.fasta --query $QUERY --outfmt 6 --out $OUT_PATH/${SMPLE}_blastp.out --block-size 18 --salltitles --quiet --more-sensitive

# test command
#diamond blastp --db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/CARD/protein_fasta_protein_homolog_model.fasta --query /xdisk/barberan/mtoucedasuarez/landuse/gene_prediction/protein_fasta/11B-T1-3_pro.fasta --outfmt 6 --out /xdisk/barberan/mtoucedasuarez/landuse/args/CARD/test_blastp.out  --block-size 0.5  --salltitles  --quiet --more-sensitive


# for loop over samples:
for sample in `awk '{print $1}' /xdisk/barberan/mtoucedasuarez/landuse/s_regex.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # translate the nucleotide fasta reads into protein sequences to be able to run diamond on them

    diamond blastx --db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/CARD/protein_fasta_protein_homolog_model.fasta --query /xdisk/barberan/mtoucedasuarez/landuse/clean/fasta/${sample}_R1_paired.fa --outfmt 6 --max-target-seqs 1 --evalue 0.000001 --out /xdisk/barberan/mtoucedasuarez/landuse/args/card/pre-assembly/${sample}_blastp.out --salltitles  --quiet --more-sensitive

done
