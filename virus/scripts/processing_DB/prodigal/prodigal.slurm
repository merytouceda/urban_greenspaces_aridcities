#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=prodigal
#SBATCH --output=prodigal.out
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=168gb
#SBATCH --time=24:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

cd /xdisk/egornish/mtoucedasuarez/vocs/

for sample in `awk '{print $1}' assembly/megahit_output/samples1.txt`
do

        if [ "$sample" == "sample" ]; then continue; fi

        ~/software/Prodigal-2.6.1/prodigal -i assembly/megahit_output/contigs/filtered/${sample}_filt_contigs.fasta -a gene_prediction/${sample}_pro.fa -d gene_prediction/${sample}_nuc.fa -f gff -p meta -q -m

done