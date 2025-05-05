#!/bin/bash -l
#SBATCH --job-name=binning
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
cd /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/mags/binning
source config.sh

# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# usage
#-a STR          metagenomic assembly file
#-o STR          output directory
#-t INT          number of threads (default=1)
#-m INT	   	 amount of RAM available (default=4)
#-l INT		 minimum contig length to bin (default=1000bp). Note: metaBAT will default to 1500bp minimum

# add -metabat2 --maxbin2 (took it out cause concoct hadn't worked and needed to repeat just with it)
metawrap binning -o $OUT_DIR/${SAMPLE}_bins -a $ASSEMBLY_DIR/${SAMPLE}_megahit_output/${SAMPLE}_contig.fa -t 48 -m 4 --concoct $READS_DIR/*.fastq
