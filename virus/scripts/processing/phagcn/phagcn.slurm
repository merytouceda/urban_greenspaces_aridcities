#!/bin/bash -l
#SBATCH --job-name=phagcn
#SBATCH --account=gornish
#SBATCH --output=output_phagcn.txt
#SBATCH --error=errors_phagcn.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu

module load anaconda
source ~/.bashrc && conda activate
conda activate phagcn

# change directory to github clone
cd /groups/barberan/software/PhaGCN_newICTV

# run phagcn
python run_Speed_up.py --contigs /xdisk/barberan/mtoucedasuarez/landuse/virus/cluster/viral_contigs_rep_seq.fasta --len 8000

#The output file is final_prediction.csv. There are three column in this csv file: "contig_name, median_file_name, prediction".
