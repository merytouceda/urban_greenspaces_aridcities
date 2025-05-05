#!/bin/bash -l
#SBATCH --job-name=mergeandfinal
#SBATCH --account=barberan
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000


# load job configuration
cd /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/mergeandfinal
source config.sh


# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $SAMPLE_LIST | tail -n 1`
# Get file names
export ASS=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# Create variable for output from inferences
MY_CH_OUT1="$VS2_OUTDIR/${ASS}_VS2_CheckV"
MY_CH_OUT2="$DVF_OUTDIR/${ASS}_DVF_CheckV"

#################### MERGING RESULTS ############################
module load R

# get coverage file 
COV="$ASSEMBLY_DIR/${SAMPLE}_cov.txt"

# Run the filtering script
cd $FINAL
if [ ! -d $FINAL ]; then
  mkdir -p $FINAL;
fi

#Rscript --no-save $SCRIPTS_DIR/config_scripts/FilterViral.R $MY_CH_OUT2/quality_summary.tsv $MY_CH_OUT1/quality_summary.tsv $COV $ASS

# stop here and see the number that comes out so that I can control the filtering depending on the number of output viral seqs

################## get final fasta ###############################

module load python/3.6

RECAP="$FINAL/${ASS}_selection1.csv"

if [ ! -d $FINAL/fasta ]; then
  mkdir -p $FINAL/fasta;
fi

python $SCRIPTS_DIR/config_scripts/get_selection_viral.py -f $RECAP -a $ASSEMBLY_DIR/contigs/$ASS -c1 $MY_CH_OUT1 -c2 $MY_CH_OUT2 -n $ASS -o $FINAL/fasta



# echo for log
echo "job done"; date
