#!/bin/bash -l
#SBATCH --job-name=vc2
#SBATCH --account=barberan
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000


module load anaconda
source ~/.bashrc && conda activate
conda activate vContact2

source $SLURM_SUBMIT_DIR/config.sh

# echo for log
echo "job started"; hostname; date

# Get sample ID
export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# Run vcontact2
vcontact2 --raw-proteins $PRODIGAL_DICT/${SAMPLE}_pro.fa --proteins-fp $G2G_DIR/${SAMPLE}_gene-to-genome.csv --db 'ProkaryoticViralRefSeq94-Merged' --output-dir $VCONTACT_DIR
