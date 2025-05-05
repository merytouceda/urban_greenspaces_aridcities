#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=coverage
#SBATCH --output=coverage.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=16:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------



# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
cd /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/coverage
source config.sh


# echo for log                                                                                                                                                                                                                                                               
echo "job started"; hostname; date


# get file name from list
export SAMPLE_ID=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# create variable of read file
PAIR1="$READS_DIR/${SAMPLE_ID}_R1_paired.fastq.gz"
PAIR2="$READS_DIR/${SAMPLE_ID}_R2_paired.fastq.gz"

#sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta

# align
$BBMAP_DIR/bbwrap.sh ref=$ASSEMBLY_DIR/${SAMPLE_ID}_megahit_output/final.contigs.fa in=$PAIR1 in2=$PAIR2 out=$ASSEMBLY_DIR/${SAMPLE_ID}_aln.sam.gz kfilter=22 subfilter=15 maxindel=80 nodisk

# coverage log
$BBMAP_DIR/pileup.sh in=$ASSEMBLY_DIR/${SAMPLE_ID}_aln.sam.gz out=$ASSEMBLY_DIR/${SAMPLE_ID}_cov.txt


# echo for log                                                                                                                                                                                                                                                               
echo "job done"; date
