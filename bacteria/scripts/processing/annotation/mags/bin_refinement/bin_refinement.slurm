#!/bin/bash -l
#SBATCH --job-name=bin_refinement
#SBATCH --account=gornish
#SBATCH --output=outputr_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu                                                                                                        
#SBATCH --partition=standard  
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=4000

# load environment
#source $CONDA/etc/profile.d/conda.sh
module load anaconda
source ~/.bashrc && conda activate
conda activate metawrap-env

# load job configuration
cd /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/mags/bin_refinement
source config.sh

# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`


# usage
#-o STR          output directory
#-t INT          number of threads (default=1)
#-m INTmemory available (default=40)
#-c INT          minimum % completion of bins [should be >50%] (default=70)
#-x INT          maximum % contamination of bins that is acceptable (default=10)

#-A STRfolder with metagenomic bins (files must have .fa or .fasta extension)
#-B STRanother folder with metagenomic bins
#-C STRanother folder with metagenomic bins

# run bin refinement
metawrap bin_refinement -o $OUT_DIR/${SAMPLE}_refined_bins -A $BINNING_DIR/${SAMPLE}_bins/metabat2_bins -B $BINNING_DIR/${SAMPLE}_bins/maxbin2_bins -C $BINNING_DIR/${SAMPLE}_bins/concoct_bins -c 70 -x 10 
