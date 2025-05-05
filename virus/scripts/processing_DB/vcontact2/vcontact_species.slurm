#!/bin/bash -l
#SBATCH --job-name=vc2
#SBATCH --account=gornish
#SBATCH --output=output_vc2_sp.txt
#SBATCH --error=errors_vc2_sp.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu

module load anaconda
source ~/.bashrc && conda activate
conda activate vContact2

# other dependencies
module load mcl/14-137
module load blast/2.10.0 
module load diamond/2.0.9

# echo for log
echo "job started"; hostname; date

#add clustal one path
export PATH=$HOME/.conda/envs/vContact2/bin:$PATH

# Run vcontact2
vcontact2 --raw-proteins /xdisk/barberan/mtoucedasuarez/urban/virus/cluster/viral_contigs_cluster_pro.fa --proteins-fp /xdisk/barberan/mtoucedasuarez/urban/virus/cluster/viral_contigs_cluster_gene-to-genome.csv --db 'ProkaryoticViralRefSeq94-Merged' --output-dir /xdisk/barberan/mtoucedasuarez/urban/virus/vcontact2/cluster_species