#!/bin/bash -l
#SBATCH --job-name=g2g
#SBATCH --account=barberan
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000

module load R

source $SLURM_SUBMIT_DIR/config.sh

# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

cd $VCONTACT_DIR
Rscript $SLURM_SUBMIT_DIR/parse_genetogenome2.R  $PRODIGAL_DIR/${SAMPLE}_nuc.fa $SAMPLE

# echo for log
echo "job done"; date
