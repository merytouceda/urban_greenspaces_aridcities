#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=parse_dvf
#SBATCH --output=parse_dvf.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00

# --------------------------------------------------                                                                                                                                          
# Load modules here                                                                                                                                                                      
# --------------------------------------------------
module load python/3.6

# --------------------------------------------------                                                                                                                                          
# Run code here                                                                                                                                                                           
# -------------------------------------------------- 
# change directory
source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/deepvirfinder/dvf_parse/config.sh
cd $DVF_OUTDIR

# Get sample ID
export ASS=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`


DVF_OUTFILE="${ASS}_gt1500bp_dvfpred.txt"

# Test on one sample
#python /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/HPC_ViralDetection/config_scripts/parse_dvf.py -f 11B-T1-3_contig.fa_gt1500bp_dvfpred.txt -a /xdisk/barberan/mtoucedasuarez/landuse/assembly/11B-T1-3_megahit_output/11B-T1-3_contig.fa -n test_megahit.fa -o parse


# run array
# parse_dvf.py usage
# -f: dvf output file
# -a: assembly output file (contigs)
# -n: name of the sample

python $PARSE_GITHUB_DIR/config_scripts/parse_dvf.py -f $DVF_OUTFILE -a $ASS_OUTDIR/$ASS -n $ASS -o $PARSE_OUTDIR
