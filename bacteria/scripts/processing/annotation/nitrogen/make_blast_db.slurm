#!/bin/bash -l
#SBATCH --job-name=blast_db
#SBATCH --account=barberan
#SBATCH --output=output_blast_db.txt
#SBATCH --error=errors_blast_db.txt
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu


# We are going to create a blast db with our viral species (output from clustering)

module load blast/2.10.0

IN_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/nitrogen/NCyc/data"

cd $IN_DIR
#gzip -d $IN_DIR/NCyc_100.faa.gz
makeblastdb -in $IN_DIR/NCyc_100.faa -out $IN_DIR/NCyc_100_db -dbtype prot -title NCyc_100 -parse_seqids -max_file_sz '4GB'
