#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=vs2
#SBATCH --output=vs2.out
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
module load anaconda
source ~/.bashrc && conda activate
conda activate vs2
# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
cd /groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/vs2 
source config.sh

#for sample in `awk '{print $1}' s_regex.txt`
#do
#   if [ "$sample" == "sample" ]; then continue; fi
   
#   virsorter run -w virus/vs2/${sample}.out -i assembly/contigs/${sample}_contig.fa --include-groups RNA -j 4 all

#done

# array of jobs: 
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

virsorter run -w $OUT_PATH/$SMPLE.out -i $IN_DIR/$SMPLE -j 4 all --min-length 1500 
