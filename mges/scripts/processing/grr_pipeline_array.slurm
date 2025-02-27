#!/bin/bash

#SBATCH --job-name=wgrr_pipeline
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=50:00:00
#SBATCH --account=barberan


module load anaconda
source ~/.bashrc && conda activate
conda activate mmseqs2
module load python/3.9/3.9.10 

# Initialize config.sh and list of files
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/urban/annotation/grr/mmseqs_method/per_sample/config.sh
export FILE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

start_time=$(date +%s)  # Capture start time
echo "Start time: $start_time"

# CLUSTER 
# Run mmseqs clustering of the protein sequences of plasmids and bacteria to remove duplicates
echo "Ok here we go, running clustering of bacteria"
mmseqs easy-cluster $BACT_PROTS_DIR/${FILE}.contigs.fa_pro.fa $BACT_PROTS_DIR/${FILE}_bact_pro_urban tmp --min-seq-id 0.95 --cov-mode 1 -c 0.7 --threads 48
echo "Done! Running clustering of plasmids"
mmseqs easy-cluster $PLAS_PROTS_DIR/${FILE}.contigs.fa_pro.fa $PLAS_PROTS_DIR${FILE}_plas_pro_urban tmp --min-seq-id 0.95 --cov-mode 1 -c 0.7 --threads 48
echo "Success!" 


# FILTER AND PREPARE
# Filter the sequencs from the bacterial protein set that have been identified as plasmids
echo "Filtering the plasmids out of the bacterial protein files"
pyhton3 $JOB_DIR/filter_fasta.py $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq.fasta $BACT_PROTS_DIR/${FILE}_bact_pro_urban_rep_seq.fasta $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta

# Rename the fasta headers so they read PLASMID_.. and CONTIG_..
echo "...changing the header names ..." 
sed 's/^>/\>CONTIG_/' $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta > $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta
sed 's/^>/\>PLASMID_/' $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq.fasta > $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta

# Concatenate the two files to create the input file for mmseqs
echo " ...and concatenating the files to get ready for alignment..."
cat $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta > $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fast




# ALL VS ALL ALIGNMENT
# Run mmseqs all against all
echo "Alright! Starting the all vs. all alginment. Good luck!" 
mmseqs easy-search $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta $MMSEQS_DIR/${FILE}_plas_all.m8 $MMSEQS_DIR/tmp -s 7.5 --threads 24 --format-output "query,target,qcov,tcov,fident,evalue,bits"  > $MMSEQS_DIR/${FILE}_plas_all.mmseqs.search.log
echo "And it worked!"



# FILTER AND WGRR CALCULATION
# Filter out: 
   # Remove same query-target
   # Keep only the best hit (smallest e-value and highest fident) per query
echo "Now I am just filtering the m8 file"
python3 $JOB_DIR/process_m8.py -f -m $MMSEQS_DIR/${FILE}_plas_bact.m8 -o $MMSEQS_DIR/${FILE}_plas_bact_filtered.txt 

# Finish filtering and calculate wGRR
   # Keep only mixed hits (plasmid-bact or bact-plasmid, but not plasmid-plasmid or bact-bact)
   # Keep only the best bi-directional hits (you are my best hit and I am yours)
   # Calculate wGRR (https://gitlab.pasteur.fr/jgugliel/wgrr)
echo " and calculating the wGRR"
python3 $JOB_DIR/calculate_grr.py -fm -fb -w -m $MMSEQS_DIR/${FILE}_plas_bact_filtered.txt -a $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta -b $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta -o $MMSEQS_DIR/${FILE}_grr.txt
echo "Success! Well done, you! Check out the output for the much deserved wGRR values!"


end_time=$(date +%s)  # Capture end time
elapsed_time=$((end_time - start_time))

echo "End time: $end_time"
echo "Processing time: $elapsed_time seconds"
