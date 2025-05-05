#!/bin/bash -l   
IN_LIST="/groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/sample_list.txt" # list of files to process  
SAMPLE_LIST="/groups/egornish/mtoucedasuarez/scripts/jobs/urban/virus/sample_id_list.txt" # I need the sample name for the coverage, this is list of samples
FINAL="/xdisk/barberan/mtoucedasuarez/urban/virus/final"
DVF_OUTDIR="/xdisk/barberan/mtoucedasuarez/urban/virus/dvf" # directory with deepvirfinder output                                                                                   
VS2_OUTDIR="/xdisk/barberan/mtoucedasuarez/urban/virus/vs2" # directory with virsorter2 output
SCRIPTS_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/HPC_ViralDetection" # path to github dir with scripts Alise made
ASSEMBLY_DIR="/xdisk/barberan/mtoucedasuarez/urban/contigs" # directory with the assembled contigs and coverage files
