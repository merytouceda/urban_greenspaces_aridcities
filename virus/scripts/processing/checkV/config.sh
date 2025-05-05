#!/bin/bash -l                                                                                                                                                                                              
# conda directory                                                                                                                                                                                           
#CONDA=""                                                                                                                                                                                                   
# INPUT LIST and  DIRECTORY                                                                                                                                                                                 
IN_LIST="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/smple_list.txt" # list of files to process                                                                                              
DVF_OUTDIR="/xdisk/barberan/mtoucedasuarez/landuse/virus/deepvirfinder" # directory with deepvirfinder output                                                                                              

VS2_OUTDIR="/xdisk/barberan/mtoucedasuarez/landuse/virus/virsorter2" # directory with virsorter2 output

CHECKVDB="/groups/egornish/mtoucedasuarez/scripts/jobs/landuse/virus/checkv-db-v1.5" # path to checkv database download 

