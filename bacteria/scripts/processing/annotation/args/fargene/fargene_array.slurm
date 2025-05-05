#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=lu_fargene_array
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=5:00:00


# --------------------------------------------------
# Load modules
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate seqkit

# --------------------------------------------------
# Execute
# --------------------------------------------------
source /groups/egornish/mtoucedasuarez/scripts/jobs/landuse/args/fargene/config.sh

export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

cd /groups/barberan/software/fargene


fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_a -o ${SAMPLE}_out/class_a_out -p 40 --force 

# class b beta-lactamases 
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_b_1_2 -o ${SAMPLE}_out/class_b_1_2_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_b_3 -o ${SAMPLE}_out/class_b_3_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# class c beta-lactamases 
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_c -o ${SAMPLE}_out/class_c_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# class d beta-lactamases   
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_d_1 -o ${SAMPLE}_out/class_d_1_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model class_d_2 -o ${SAMPLE}_out/class_d_2_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# qnr
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model qnr -o ${SAMPLE}_out/qnr_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# tetracycline
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model tet_efflux -o ${SAMPLE}_out/tet_efflux_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model tet_rpg -o ${SAMPLE}_out/tet_rpg_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model tet_enzyme -o ${SAMPLE}_out/tet_enzyme_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# erm 23S RNA
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/erm_type_a.hmm -o ${SAMPLE}_out/erm_type_a_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/erm_type_f.hmm -o ${SAMPLE}_out/erm_type_f_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# mph
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/mph.hmm -o ${SAMPLE}_out/mph_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

#AAC (acetyltranferases)
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_a.hmm -o ${SAMPLE}_out/aminoglycoside_model_a_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_b.hmm -o ${SAMPLE}_out/aminoglycoside_model_b_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_c.hmm -o ${SAMPLE}_out/aminoglycoside_model_c_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_d.hmm -o ${SAMPLE}_out/aminoglycoside_model_d_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_e.hmm -o ${SAMPLE}_out/aminoglycoside_model_e_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_f.hmm -o ${SAMPLE}_out/aminoglycoside_model_f_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

# APH (aminoglycoside phosphotransferases)
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_g.hmm -o ${SAMPLE}_out/aminoglycoside_model_g_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_h.hmm -o ${SAMPLE}_out/aminoglycoside_model_h_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force
fargene --infiles $READS_DIR/${SAMPLE}_*.fastq --meta --hmm-model /groups/barberan/software/fargene/fargene_analysis/models/aminoglycoside_model_i.hmm -o ${SAMPLE}_out/aminoglycoside_model_i_out -p 40 --rerun --fasta-dir ${SAMPLE}_out/class_a_out/tmpdir/ --force

