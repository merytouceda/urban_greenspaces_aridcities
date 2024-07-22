# Soil metagenomics of urban greenspaces in arid cities
Repository for the publication: "Differences in the genomic potential of soil bacterial and phage communities between urban greenspaces and natural arid soils". 


#### This repository
This repository harbours the processing and analysis information, and final data products for the study of soil bacteria and viruses in urban greenspaces and natural soils around the arid city of Tucson (Arizona, USA). 
- Bacteria
    - gene centric: assembly, gene prediction and annotation with :
        - KEGG
        - NCyc (nitrogen cycling genes) 
        - CAZY (carbon cycling genes) 
        - CARD (anbtibiotic resistance genes)
        - BACmet (heavy metal resistance genes).
    - taxonomy: taxonomic classification of reads using Kraken2. 
    - 
    - args: annotation of antibiotic resistence genes 

- Virus
    - inference
    - annotation 




#### Directory Structure: 
-bacteria
    - genecentric 
        - README.md
        - scripts
        - data
        - pipeline.md
        - analyses
    - taxonomy
        - scripts
        - analyses
    - analyses
- virus 
    - README.md
    - inference
        - deepvirfinder
        - virsorter2
    - annotation
        - phagcn
        - phatyp
        - genomad
        - iphop
    - analyses

  


#### Field sampling design
24 samples
2 land uses, 4 vegetation types, 3 replicates per site:
   
| landuse        | site              | Vegetation      |
| -------------  | ----------------  |---------------- |
| urban          | Reid Park (RP)    | Grassland       |
| urban          | Himmel Park (HP)  | Grassland       |
| urban          | Rose Canyon (RC)  | Forest          |
| natural        | Sabino Canyon (SC)| Shrubland       |
| natural        | SRER-11B          | Grassland       |
| natural        | SRER-8            | Grassland       |
| natural        | SRER-Ex45         | Grassland       |
| natural        | SRER-UAB          | Grassland       |

 *SRER stands for Santa Rita Experimental Range
 Replicates are taken in a transect at 3,6, and 9 meters for each site.
   
   
