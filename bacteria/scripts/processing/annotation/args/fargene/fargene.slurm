#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=lu_fargene_test
#SBATCH --output=lu_fargene_test.out
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=10:00:00


# --------------------------------------------------
# Load modules
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate seqkit

# --------------------------------------------------
# Execute
# --------------------------------------------------

# test on on sample only

TEST_IN_DIR="/xdisk/barberan/mtoucedasuarez/landuse/clean/uncompressed/test"
TEST_OUT_DIR="/xdisk/barberan/mtoucedasuarez/landuse/args/fargene"

cd /groups/barberan/software/fargene

#class a beta-lactamases 
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_a -o $TEST_OUT_DIR/class_a_out -p 40 --force 

# class b beta-lactamases 
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_b_1_2 -o $TEST_OUT_DIR/class_b_1_2_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_b_3 -o $TEST_OUT_DIR/class_b_3_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# class c beta-lactamases 
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_c -o $TEST_OUT_DIR/class_c_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# class d beta-lactamases   
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_d_1 -o $TEST_OUT_DIR/class_d_1_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model class_d_2 -o $TEST_OUT_DIR/class_d_2_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# qnr
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model qnr -o $TEST_OUT_DIR/qnr_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# tetracycline
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model tet_efflux -o $TEST_OUT_DIR/tet_efflux_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model tet_rpg -o $TEST_OUT_DIR/tet_rpg_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model tet_enzyme -o $TEST_OUT_DIR/tet_enzyme_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# erm 23S RNA
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model erm_type_a -o $TEST_OUT_DIR/erm_type_a_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model erm_type_f -o $TEST_OUT_DIR/erm_type_f_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# mph
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model mph -o $TEST_OUT_DIR/mph_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

#AAC (acetyltranferases)
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_a -o $TEST_OUT_DIR/aminoglycoside_model_a_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_b -o $TEST_OUT_DIR/aminoglycoside_model_b_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_c -o $TEST_OUT_DIR/aminoglycoside_model_c_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_d -o $TEST_OUT_DIR/aminoglycoside_model_d_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_e -o $TEST_OUT_DIR/aminoglycoside_model_e_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_f -o $TEST_OUT_DIR/aminoglycoside_model_f_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

# APH (aminoglycoside phosphotransferases)
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_g -o $TEST_OUT_DIR/aminoglycoside_model_g_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_h -o $TEST_OUT_DIR/aminoglycoside_model_h_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force
fargene --infiles $TEST_IN_DIR/*.fastq --meta --hmm-model aminoglycoside_model_i -o $TEST_OUT_DIR/aminoglycoside_model_i_out -p 40 --rerun --fasta-dir $TEST_OUT_DIR/class_a_out/tmpdir/ --force

