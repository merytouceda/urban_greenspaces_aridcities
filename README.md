# LandUse
Repository for the land use project. Understanding the effects of land use practices on arid soil microbial communities. 


#### Experimental design
Shotgun metagenomics data
24 samples
4 land uses, 2 sites per land use, 3 replicates per site:
   
   **land use**     **site1**             **site 2**
   urban            Reid Park (RP)        Himmel Park (HP)
   natural          Sabino Canyon (SC)    Rose Canyon (RC)
   old pasture      SRER*-11B             SRER-8
   recent pasture   SRER-UAB              SRER-Ex45
 
 *SRER stands for Santa Rita Experimental Range
 Replicates are taken in a transect at 3,6, and 9 meters. 
   
   
   
#### Contents
This repository holds the different processeing and subsequent analyses performed on this dataset: 
- gene centric: assembly, gene prediction and annotation.  
- virus: assembly of reads and viral sequence inference from contigs. 
- taxonomy: taxonomic classification of reads using Kraken2. 
- args: annotation of antibiotic resistence genes, 



#### Directory Structure: 
- genecentric 
    - scripts
    - data
    - README.md
    - pipeline.md
    - analyses
- virus 
    - deepvirfinder
    - virsorter2
    - README.md
    - analyses
- taxonomy
    - scripts
    - analyses
- args
    - scripts
    - analyses
  
