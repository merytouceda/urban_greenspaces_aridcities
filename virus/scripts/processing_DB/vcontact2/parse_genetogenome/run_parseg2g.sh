#!/bin/bash -l

# load job configuration
#IN_LIST="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/smple_list.txt"
source /groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/vcontact2/parse_genetogenome/config.sh

#
#makes sure sample file is in the right place
#
if [[ ! -f "$IN_LIST" ]]; then
    echo "$IN_LIST does not exist. Please provide the path for a list of datasets to process. Job terminated."
    exit 1
fi

# get number of samples to process
export NUM_JOB=$(wc -l < "$IN_LIST")

# submit co_assemblies
echo "launching parse_genetogenome.slurm"

JOB_ID=`sbatch --job-name g2g -a 1-$NUM_JOB parse_genetogenome.slurm`
