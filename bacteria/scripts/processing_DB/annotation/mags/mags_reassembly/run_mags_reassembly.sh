#!/bin/bash -l

# load job configuration
#IN_LIST="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/smple_list.txt"
source /groups/egornish/mtoucedasuarez/scripts/jobs/urban/mags/mags_reassembly/config.sh

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
echo "launching mags_reassembly.slurm as a job."

JOB_ID=`sbatch --job-name mags_reassembly -a 1-$NUM_JOB mags_reassembly.slurm`
