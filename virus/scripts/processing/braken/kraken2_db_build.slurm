#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=kraken2
#SBATCH --output=kraken2_db.out
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=10:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate kraken2



# --------------------------------------------------
# Run command here
# --------------------------------------------------

#kraken2-build --standard --threads 24 --db /xdisk/egornish/mtoucedasuarez/databases/kraken2/kraken2 

/home/u14/mtoucedasuarez/software/kraken2/kraken-build --download-library bacteria --db /xdisk/egornish/mtoucedasuarez/databases/kraken2/bacteria
