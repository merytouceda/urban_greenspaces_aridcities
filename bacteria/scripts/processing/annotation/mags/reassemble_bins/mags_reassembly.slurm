#!/bin/bash -l
#SBATCH --job-name=mags_reassembly
#SBATCH --account=gornish
#SBATCH --output=outputr_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu                                                                                                        
#SBATCH --partition=standard  
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000


# load job configuration
cd /groups/egornish/mtoucedasuarez/scripts/jobs/urban/mags/mags_reassembly
source config.sh

# load environment
#source $CONDA/etc/profile.d/conda.sh
module load anaconda
source ~/.bashrc && conda activate
conda activate metawrap-env


# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# usage
#-b STR          folder with metagenomic bins
#-o STR          output directory
#-1 STR          forward reads to use for reassembly
#-2 STR          reverse reads to use for reassembly

metawrap reassemble_bins -o $REASSEMBLY_DIR/${SAMPLE}_reassembly -1 $READS_DIR/${SAMPLE}_clean_1.fastq.gz -2 $READS_DIR/${SAMPLE}_clean_2.fastq.gz -m 35000 -t 10 -c 70 -x 10 -b $BINS_DIR
