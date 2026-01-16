######################################################################################################
# average genome size(AGS) and average 16S rRNA copy number (ACN)
######################################################################################################

# output should be exported to a trait folder
# output of AGS should be exported to trait/ags/
# output of ACN should be exported to trait/acn/
mkdir trait
mkdir -p trait/ags/
mkdir -p trait/acn/

ags.sh

module load anaconda
source /groups/jneilson/chenyj/conda/etc/profile.d/conda.sh

conda activate sortmerna4.3.1

# Note1: The following will generate single_cogs_count.tsv and total base pairs,
# then R is used to calculate the final AGS result (see below)
# Note2: FGSpp will crash occasionally for some samples
# Note3: The following will also generate blast output of sortmerna,
# then customized command is used to parse the blast output (see below)

for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # convert R1 and R2 to interleave files
    ~/software/bbmap/bbduk.sh \
    in1=/your/path/${sample}_clean_1.fastq.gz \
    in2=/your/path/${sample}_clean_2.fastq.gz \
    out=/your/path/${sample}_interleave.fa \

    # rename R1 and R2
    ~/software/bbmap/reformat.sh \
    in1=/your/path/${sample}_1.fa \
    in2=/your/path/${sample}_2.fa \
    out=/your/path/${sample}_interleave.fa \
    trimreaddescription=t \
    addslash=t \
    spaceslash=f

    # run sortmerna
    sortmerna \
    --ref /xdisk/jneilson/chenyj/smr_v4.3_default_db.fasta \
    --reads /your/path/${sample}_interleave.fa \
    --workdir /your/path/${sample}_tmp \
    --fastx \
    -m 3072 \
    -e 1e-1 \
    --blast 1 \
    --num_alignments 1 \
    -threads 94

    # copy blast file to acn
    mv /your/path/${sample}_tmp/out/aligned.blast /your/path/${sample}.blast

    # remove intermediate files
    rm -R /your/path/${sample}_tmp
    rm /your/path/${sample}_1.fa /your/path/${sample}_2.fa

    # Average genome size
    # First, use fragescanplusplus to predict orf on raw reads, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh
    # you might want to change the training file, i.e., illumina_5 if FGSpp crashes
    ~/software/FragGeneScanPlusPlus/FGSpp \
    -s /your/path/trait/${sample}_interleave.fa \
    -o /your/path/trait/ags/${sample}_fgs \
    -w 0 \
    -r ~/software/FragGeneScanPlusPlus/train/ \
    -t illumina_5 \
    -p 94 \
    -c 400

    # Second, calculate the total number of base pairs
    egrep -v "^>" /your/path/trait/${sample}_interleave.fa | wc | awk '{print $3-$1}' > /your/path/trait/ags/${sample}_totalbp.txt

    # Third, annotate single copy genes using uproc, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh
    ~/software/uproc/bin/uproc-prot \
    --threads 94 \
    --output /your/path/trait/ags/${sample}_uproc \
    --preds \
    --pthresh 3 \
    /groups/barberan/metagenome_trait/ags_acn/uproc_scg_db/SINGLE_COPY_COGS_DB \
    /groups/barberan/metagenome_trait/ags_acn/uproc_scg_db/modeldir/model \
    /your/path/trait/ags/${sample}_fgs.faa

    # Fourth, get single copy gene count, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh
    # This will generate ${sample}_single_cogs_count.tsv, we need this and total number of base pairs to calculate the AGS (see below)
    cut -f1,3,4,5 -d"," /your/path/trait/ags/${sample}_uproc | \
    awk 'BEGIN {OFS="\t"; FS=","} {
    if (array_score[$1]) {
    if ($4 > array_score[$1]) {
      array_score[$1] = $4
      array_line[$1, $3] = $2
    }
    } else {
    array_score[$1] = $4
    array_line[$1, $3] = $2
    }
    } END {
    printf "%s\t%s\n", "cog","cov"
    for (combined in array_line) {
    split(combined, separate, SUBSEP)
    array_length[separate[2]]= array_length[separate[2]] + array_line[combined]
    }
    for ( c in array_length ) {
    printf "%s\t%s\n", c,array_length[c]
    }
    }' > /your/path/trait/ags/${sample}_single_cogs_count.tsv

    # remove intermediate files
    rm /your/path/${sample}_interleave.fa
done # change me to your path

# AGS estimation. From above analyses, we already have the sample_single_cogs_count.tsv and the sample_totalbp.txt
# single_cogs_count.tsv can be used to estimate the number of genomes, then AGS is
# calculated as the number of genomes divided by the sample_totalbp.txt
# All codes are copied from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/

ags.R

# ACN estimation. From the above analyses, we already have the sample.blast
# Now we need to calculate the 16S rRNA coverage using sample.blast, then
# ACN is calculated as 16S rRNA coverage divided by the number of genomes
# all codes are copied from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/acn.sh

acn_sortmerna.R

for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # filter blast output from sortmerna, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/acn.sh
    # l: minimum length to be used in 16S rRNA filtering (default 30)
    # i: minimum identity to be used in 16S rRNA filtering (default 85)
    # e: e-value to be used in 16S rRNA filtering (default 1e-15)
    # s: 16S rRNA reference length (default E. cloi: 1542bp)

    awk -v l=30 -v i=85 \
    -v e=1e-15 -v s=1542 \
      '{
        if ( $11 <= e && $4 >= l && $3 >= i ) {
        n_nuc = $10 -$9 +1;
        n_nuc = sqrt(n_nuc*n_nuc)
        n_nuc_tot = n_nuc + n_nuc_tot
      }
    } END {
      print n_nuc_tot/s
    }' trait/acn/${sample}.bac_archaea > trait/acn/${sample}_16S_coverage.txt

done

acn.R



######################################################################################################
# GC content and GC variance
######################################################################################################

# we calculate the gc content and gc variance of raw reads

# output of should be exported to trait/gc
mkdir -p trait/gc

module load anaconda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate seqkit

for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # convert R1 and R2 to interleave files
    ~/software/bbmap/bbduk.sh \
    in1=/your/path/${sample}_clean_1.fastq.gz \
    in2=/your/path/${sample}_clean_2.fastq.gz \
    out=trait/gc/${sample}_interleave.fa \

    seqkit fx2tab \
    /your/path/${sample}_interleave.fa \
    -g -i -n -j 94 > trait/gc/${sample}_gc.txt

    # mean GC
    awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' trait/gc/${sample}_gc.txt > trait/gc/${sample}_gc_mean.txt

    # GC standard deviation
    awk '{sum+=$2; sumsq+=$2*$2} END {print sqrt(sumsq/NR - (sum/NR)**2)}' trait/gc/${sample}_gc.txt > trait/gc/${sample}_gc_var.txt

    rm trait/gc/${sample}_gc.txt
    rm trait/gc/${sample}_interleave.fa

done





######################################################################################################
# Growth rate
######################################################################################################

# gRodon is used to calculate growth rate: https://github.com/jlw-ecoevo/gRodon

# Growth rate prediction requires ribosomal protein genes
# The first step is to annotate ribosomal protein genes in each sample using kofamscan
# See here for a list of KOs associated with ribosomal protein gene: https://www.genome.jp/kegg-bin/show_brite?ko03011.keg

growth_kofamscan.pbs

module load anaconda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate kofamscan

for sample in `awk '{print $1}' samples1.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    /groups/barberan/metagenome_annotation/KEGG_kofamscan/kofamscan/bin/kofam_scan-1.3.0/exec_annotation -f mapper \
    -p /groups/barberan/metagenome_annotation/KEGG_kofamscan/profiles/rb.hal \
    -k /groups/barberan/metagenome_annotation/KEGG_kofamscan/ko_list \
    --cpu 28 \
    --tmp-dir trait/growth/rb_annotation/tmp \
    -o trait/growth/rb_annotation/${sample}_rb_kofamscan_raw.txt \
    nucleotide_protein/${sample}_pro.fa

    echo -e "gene\tKO" | cat - trait/growth/rb_annotation/${sample}_rb_kofamscan_raw.txt > trait/growth/rb_annotation/${sample}_rb_kofamscan.txt
    rm trait/growth/rb_annotation/${sample}_rb_kofamscan_raw.txt
done


# calculate gene coverage, gene coverage can be taken into account when calculating
# growth rate, see details in https://github.com/jlw-ecoevo/gRodon
growth_coverm.pbs

module load anaconda
source ~/miniconda3/etc/profile.d/conda.sh

conda activate coverm

for sample in `awk '{print $1}' samples1.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    TMPDIR=. coverm contig \
    -1 /your/path/clean_reads/${sample}_clean_1.fastq \
    -2 /your/path/clean_reads/${sample}_clean_2.fastq \
    -r /your/path/${sample}_nuc.fa \
    -p bwa-mem \
    -t 28 \
    -m mean \
    --min-read-percent-identity 0.95 \
    --min-read-aligned-length 45 \
    -o trait/growth/gene_coverage/${sample}_gene_coverage.txt
done # change me to your path

# Everthing is ready to go, then the gRodon pakcage is used to calculate growth rate
# This program needs a set of genes (i.e., prodigal output) and a set of ribosomal protein
# genes (i.e., the output of the above hmmsearch and kofamscan). It calculates the community-level
# growth rate of each sample.

module load R

R --file=growth_pred.R
