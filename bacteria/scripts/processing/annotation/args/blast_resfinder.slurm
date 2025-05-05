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
#source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/config.sh

#export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
#$QUERY=$IN_DIR/$SMPLE
# array command
#blast -d $DB_DIR/protein_fasta_protein_homolog_model.fasta -q $QUERY --outfmt 6 --out $OUT_PATH/${SMPLE}_blastp.out --block-size 18 --salltitles --quiet --more-sensitive

# test command
#makeblastdb -in /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/resfinder_db/all.fsa -blastdb_version 5 -title "resfinder" -dbtype nucl

#blastn -db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/resfinder_db/all.fsa -query /xdisk/barberan/mtoucedasuarez/landuse/gene_prediction/11B-T1-3_prodigal_output/11B-T1-3_gene.fasta -outfmt 6 -out /xdisk/barberan/mtoucedasuarez/landuse/args/resfinder/test_blast.out   


# for loop over samples:
for sample in `awk '{print $1}' /xdisk/barberan/mtoucedasuarez/landuse/s_regex.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    # convert reads from fastq to fasta to be able to blast
    #/home/u14/mtoucedasuarez/software/seqtk/seqtk seq -a /xdisk/barberan/mtoucedasuarez/landuse/clean/${sample}_R1_paired.fastq.gz > /xdisk/barberan/mtoucedasuarez/landuse/clean/fasta/${sample}_R1_paired.fa

    blastn -db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/resfinder_db/all.fsa -query /xdisk/barberan/mtoucedasuarez/landuse/clean/fasta/${sample}_R1_paired.fa -outfmt 6 -max_target_seqs 1 -evalue 0.000001 -out /xdisk/barberan/mtoucedasuarez/landuse/args/resfinder/pre-assembly/${sample}_blast.out #try to thread it

done


#blastn -db /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/databases/resfinder_db/all.fsa -query /xdisk/barberan/mtoucedasuarez/landuse/gene_prediction/11B-T1-3_gene.fasta -outfmt 6 -out /xdisk/barberan/mtoucedasuarez/landuse/args/resfinder/11-T1-3_blast.out
