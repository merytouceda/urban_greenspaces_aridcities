#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_mags
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load samtools
# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

source $SLURM_SUBMIT_DIR/config.sh

# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

    # mapping
~/software/bwa/bwa mem $MAPPING_DIR/mags_index $READS_DIR/${SAMPLE}_R1_paired.fastq.gz $READS_DIR/${SAMPLE}_R2_paired.fastq.gz -t 28 > $MAPPING_DIR/${SAMPLE}.sam
    
# convert sam to bam
samtools view -@ 28 -F 4 -bo $MAPPING_DIR/${SAMPLE}_raw.bam $MAPPING_DIR/${SAMPLE}.sam
# sort bam file
samtools sort -@ 28 -m 6G -o $MAPPING_DIR/${SAMPLE}.bam $MAPPING_DIR/${SAMPLE}_raw.bam
# index bam file, this is optional
samtools index $MAPPING_DIR/${SAMPLE}.bam

# remove intermediate files
rm $MAPPING_DIR/${SAMPLE}.sam
rm $MAPPING_DIR/${SAMPLE}_raw.bam



# echo for log
echo "job done"; date