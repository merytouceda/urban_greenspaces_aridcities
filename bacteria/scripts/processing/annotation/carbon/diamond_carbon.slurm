#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=diamondcazy
#SBATCH --output=diamondcazy.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=20:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load diamond/2.0.9

# for loop over samples:
for sample in `awk '{print $1}' /xdisk/barberan/mtoucedasuarez/landuse/s_regex.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    diamond blastx --db /xdisk/egornish/mtoucedasuarez/databases/CAZyDB.09242021.fa --query /xdisk/barberan/mtoucedasuarez/landuse/clean/fasta/${sample}_R1_paired.fa --outfmt 6 --out /xdisk/barberan/mtoucedasuarez/landuse/carbon/${sample}_blastx.out  --max-target-seqs 1 --evalue 0.000001 --salltitles  --quiet --more-sensitive

done
