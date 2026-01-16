

This directory contains scripts and workflows for analyzing mobile genetic elements in urban greenspaces and natural arid soils.

## Directory structure
mges/
├── data/                    # Data files for MGE analysis 
├── scripts/
│   ├── processing/         # Data processing and filtering scripts
│   └── analyses/           # Statistical analyses and visualization scripts
└── README.md


## Scripts
### Processing Scripts

Located in `scripts/processing/`, these scripts are used for a series of operations needed to obtain final analysis products

- **`grr_pipeline_array.slurm`** - SLURM array job script for running the Gene Relatedness Ratio (GRR) pipeline at scale
- **`calculate_grr.py`** - Calculates gene relatedness ratios
- **`retrieve_mge_contigs.py`** - Extracts contigs identified as mobile genetic elements
- **`retrieve_mge_genes.py`** - Retrieves genes associated with identified MGEs
- **`retrieve_plasmid_contigs.py`** - Extracts plasmid contigs from assemblies
- **`process_m8.py`** - Processes BLAST m8 format output
- **`filter_count_table.py`** - Filters count tables for downstream analysis
- **`filter_fasta.py`** - Filters FASTA sequences based on specified criteria
- **`fasta_head_tail_retrieve.py`** - Retrieves header and tail sequences from FASTA files


## Analysis Scripts

Located in `scripts/analyses/`, these R scripts perform statistical analyses and generate visualizations.

### MGE Characterization

- **`mge_types_plot.R`** - Visualizes the distribution and types of MGEs across samples
- **`wGRR_viz.R`** - Visualizes weighted Gene Relatedness Ratio (wGRR) results
- **`host_range.R`** - Analyzes and visualizes the host range of identified MGEs
- **`host_taxonomy.R`** - Examines taxonomic composition of MGE host organisms

### Functional Analysis

- **`amrfinder_abricate.R`** - Analyzes antimicrobial resistance genes found in plasmids using AMRFinder and Abricate results
- **`bacmet_ptus.R`** - Analyzes BacMet (metal resistance) found in plasmids
- **`ptus_analyses.R`** - Ecological analysis of plasmid taxonomic units (PTUs), or plasmid sequences
- **`kegg_plasmids_analysis.R`** - KEGG functional annotation analysis for plasmid-associated genes
- **`unk_orf_counts.R`** - Quantifies and analyzes unknown ORFs (open reading frames)

### Defense Systems

- **`crispr_investment.R`** - Analyzes CRISPR-Cas systems and bacterial defense investment strategies

## Data Directory

The `data/` contains:
- Processed plasmid contigs and annotations
- Count tables
- Metadata files linking samples to environmental categories