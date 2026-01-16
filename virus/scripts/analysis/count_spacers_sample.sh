#!/bin/bash

# Output file
output_file="crispr_spacer_counts.tsv"

# Header for the output file
echo -e "Sample\tNum_Spacers" > "$output_file"

# Loop through all fasta files
for file in *_crispr_spacers.fa; do
    # Extract sample name from filename
    sample_name=$(basename "$file" _crispr_spacers.fa)

    # Count the number of sequences (lines starting with '>')
    seq_count=$(grep -c "^>" "$file")

    # Append to output file
    echo -e "$sample_name\t$seq_count" >> "$output_file"
done

echo "Results saved in $output_file"
