#!/bin/bash -l
#SBATCH --job-name=phatype
#SBATCH --account=gornish
#SBATCH --output=output_phatype.txt
#SBATCH --error=errors_phatype.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate phatyp
# --------------------------------------------------
# Execute commands here
# --------------------------------------------------
cd /groups/barberan/software/PhaTYP 

MID_DIR="/xdisk/barberan/mtoucedasuarez/landuse/virus/phatype/tmp"
OUT_DIR="/xdisk/barberan/mtoucedasuarez/landuse/virus/phatype"
IN_FILE="/xdisk/barberan/mtoucedasuarez/landuse/virus/cluster/viral_contigs_rep_seq.fasta"

python preprocessing.py --contigs $IN_FILE --len 1500 --midfolder $MID_DIR
python PhaTYP.py --out $OUT_DIR/landuse_phatype_out_1500.csv --midfolder $MID_DIR
